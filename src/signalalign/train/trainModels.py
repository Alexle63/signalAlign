#!/usr/bin/env python
"""Train HMMs for alignment of signal data from the MinION
"""

import sys
import os
import urllib.parse
import textwrap
import yaml
import h5py

from argparse import ArgumentParser
from random import shuffle
from shutil import copyfile

from multiprocessing import Process, current_process, Manager

from signalalign import parseFofn
from signalalign.signalAlignment import multithread_signal_alignment_samples, \
    create_signalAlignment_args, SignalAlignSample
from signalalign.hiddenMarkovModel import get_model
from signalalign.utils import processReferenceFasta
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.bwaWrapper import buildBwaIndex
from py3helpers.utils import merge_lists


def generateConfig(config_path):
    if os.path.exists(config_path):
        raise RuntimeError
    config_content = textwrap.dedent("""\
                # SignalAlign model training config file
                output_dir: signalAlign_unittest/
                samples: [{
                    fast5_dir: ../tests/minion_test_reads/canonical_ecoli_R9/,
                    fofn:,
                    positions_file:,
                    motif:,
                    label:,
                }]
                bwa_reference: ../tests/test_sequences/E.coli_k12.fasta
                bwt:
                stateMachineType: threeState
                in_T_Hmm: ../models/testModelR9_5mer_acegot_template.model
                in_C_Hmm: ../models/testModelR9_5mer_acegot_complement.model
                templateHdp:
                complementHdp:
                iterations: 3
                training_bases: 10000
                job_count: 4
                diagonal_expansion:
                constraint_trim:
                twoD: true
                alignment_file:
                DEBUG:
                test: true
                """)
    fH = open(config_path, "w")
    fH.write(config_content)
    fH.flush()
    fH.close()


def trainHMM(samples, working_folder, config, transitions=True, emissions=False):
    """Train HMM transitions.

    :param samples: list of SignalAlignSamples
    :param config: config dictionary with required arguments
    :param working_folder: FolderHandler with directory created
    :param transitions: boolean option to train transitions
    :param emissions: boolean option to train normal emissions

    :return: trained HMM model files [template_hmm, complement_hmm, template_hdp, complement_hdp]

    config required keys
    _____________
    output_dir : path to output directory
    samples: list of dictionary of samples
            fast5_dir:  path to fast5s,
            fofn: file with paths to fast5s
            positions_file: changes to nucleotides by position,
            edited_fw_reference: forward reference if different than canonical
            edited_bw_reference: backward reference if needed

    path_to_bin: path to signalMachine
    in_T_Hmm: path to template HMM model file
    in_C_Hmm: path to complement HMM model file
    templateHdp: path to template HDP model file
    complementHdp: path to complement HDP model file
    twoD: boolean option for 2D reads
    training_bases: number of bases to use for each sample during training
    job_count: number of processes to use when running SignalAlign
    iterations: number of iterations over the entire updating pipeline
    stateMachineType: type of stateMachine ("threeStateHdp" or "threeStateHmm")
    diagonal_expansion: alignment algorithm param to expand how far a path can get from guide alignment
    constraint_trim: alignment algorithm param for how much to trim the guide alignment anchors
    """
    if transitions is False:
        assert emissions is True, "Must select either transitions or emissions to train"

    # set all variables from config file
    alignment_file = config["alignment_file"]
    training_amount = config["training_bases"]
    workers = config["job_count"]
    iterations = config["iterations"]
    twoD = config['twoD']
    state_machine_type = config["stateMachineType"]
    template_model_path = config["in_T_Hmm"]
    complement_model_path = config["in_C_Hmm"]
    original_template_hdp_path = config["templateHdp"]
    original_complement_hdp_path = config["complementHdp"]
    diagonal_expansion = config["diagonal_expansion"]
    constraint_trim = config["constraint_trim"]
    test = config["test"]
    path_to_bin = config["path_to_bin"]

    # get working folder path and check if there is an alignment file
    working_folder_path = working_folder.path
    if alignment_file:
        alignment_file = os.path.abspath(alignment_file)

    # find model files
    # make some paths to files to hold the HMMs
    complement_model = None
    complement_hmm = None

    if twoD:
        assert os.path.exists(complement_model_path), \
            "Missing complement model %s" % (complement_model_path)
        complement_model_path = os.path.abspath(complement_model_path)
        complement_model = get_model(state_machine_type, complement_model_path)
        complement_hmm = working_folder.add_file_path("complement_trained.hmm")
        copyfile(complement_model_path, complement_hmm)
        assert os.path.exists(complement_hmm), "Problem copying default model to {}".format(complement_hmm)

    assert os.path.exists(template_model_path), \
        "Missing template model %s" % (template_model_path)
    template_model_path = os.path.abspath(template_model_path)
    template_model = get_model(state_machine_type, template_model_path)
    template_hmm = working_folder.add_file_path("template_trained.hmm")
    copyfile(os.path.abspath(template_model_path), template_hmm)
    assert os.path.exists(template_hmm), "Problem copying default model to {}".format(template_hmm)
    # determine if we train HMM emissions or HDP
    update_template_hmm_emissions = False
    update_complement_hmm_emissions = False
    # get the input HDP, if we're using it
    if state_machine_type == "threeStateHdp":
        assert os.path.exists(original_template_hdp_path), "Templace HDP path not found {}".format(
            original_template_hdp_path)
        original_template_hdp_path = os.path.abspath(original_template_hdp_path)
        template_hdp = working_folder.add_file_path("%s" % os.path.basename(original_template_hdp_path))
        copyfile(original_template_hdp_path, template_hdp)
        if twoD:
            assert os.path.exists(original_complement_hdp_path), "Templace HDP path not found {}".format(
                original_complement_hdp_path)
            original_complement_hdp_path = os.path.abspath(original_complement_hdp_path)

            complement_hdp = working_folder.add_file_path("%s" % os.path.basename(original_complement_hdp_path))
            copyfile(original_complement_hdp_path, complement_hdp)
        else:
            update_complement_hmm_emissions = True
            complement_hdp = None
    else:
        if emissions:
            update_template_hmm_emissions = True
            update_complement_hmm_emissions = True
        template_hdp = None
        complement_hdp = None

    # collect signal alignment arguments
    alignment_args = create_signalAlignment_args(
        destination=working_folder_path,
        stateMachineType=state_machine_type,
        in_templateHmm=template_hmm,
        in_complementHmm=complement_hmm,
        in_templateHdp=template_hdp,
        in_complementHdp=complement_hdp,
        diagonal_expansion=diagonal_expansion,
        constraint_trim=constraint_trim,
        twoD_chemistry=twoD,
        alignment_file=alignment_file,
        get_expectations=True,
        path_to_bin=path_to_bin)

    i = 0
    # start iterating
    while i < iterations:
        # align all the samples
        samples = multithread_signal_alignment_samples(samples, alignment_args, workers, trim=training_amount)
        all_sample_files = merge_lists([sample.analysis_files for sample in samples])
        assert len(all_sample_files) > 0, "Something failed in multithread signal alignment. We got no sample files"
        # load then normalize the expectations
        template_expectations_files = [x for x in all_sample_files
                                       if x.endswith(".template.expectations.tsv")]

        if len(template_expectations_files) > 0:
            template_model.add_and_normalize_expectations(files=template_expectations_files,
                                                          hmm_file=template_hmm,
                                                          update_transitions=True,
                                                          update_emissions=update_template_hmm_emissions)
        if twoD:
            complement_expectations_files = [x for x in all_sample_files
                                             if x.endswith(".complement.expectations.tsv")]
            if len(complement_expectations_files) > 0:
                complement_model.add_and_normalize_expectations(files=complement_expectations_files,
                                                                hmm_file=complement_hmm,
                                                                update_transitions=True,
                                                                update_emissions=update_complement_hmm_emissions)

        # log the running likelihood
        if len(template_model.running_likelihoods) > 0 and \
                (twoD and len(complement_model.running_likelihoods)) > 0:
            print("[trainModels_HMM] {i}| {t_likelihood}\t{c_likelihood}".format(
                t_likelihood=template_model.running_likelihoods[-1],
                c_likelihood=complement_model.running_likelihoods[-1],
                i=i))
            if test and (len(template_model.running_likelihoods) >= 2) and \
                    (config["twoD"] and len(complement_model.running_likelihoods) >= 2):
                assert (template_model.running_likelihoods[-2] < template_model.running_likelihoods[-1]) and \
                       (complement_model.running_likelihoods[-2] < complement_model.running_likelihoods[-1]), \
                    "Testing: Likelihood error, went up"
        elif len(template_model.running_likelihoods) > 0:
            print("[trainModels_HMM] {i}| {t_likelihood}".format(t_likelihood=template_model.running_likelihoods[-1],
                                                                 i=i))
            if test and (len(template_model.running_likelihoods) >= 2):
                assert (template_model.running_likelihoods[-2] < template_model.running_likelihoods[
                    -1]), "Testing: Likelihood error, went up"

        i += 1

    # if we're using HDP, trim the final Hmm (remove assignments)

    print("trainModels - finished training routine", file=sys.stdout)
    print("trainModels - finished training routine", file=sys.stderr)

    return [template_hmm, complement_hmm, template_hdp, complement_hdp]


def main():
    def parse_args():
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest="command")

        # parsers for running the full pipeline
        run_parser = subparsers.add_parser("run", help="runs full workflow ")
        run_parser.add_argument('--config', default='trainModels-config.yaml', type=str,
                                help='Path to the (filled in) config file, generated with "generate".')
        subparsers.add_parser("generate", help="generates a config file for your run, do this first")
        return parser.parse_args()

    args = parse_args()
    if args.command == "generate":
        try:
            config_path = os.path.join(os.getcwd(), "trainModels-config.yaml")
            generateConfig(config_path)
        except RuntimeError:
            print("Using existing config file {}".format(config_path))
            pass
    elif args.command == "run":
        if not os.path.exists(args.config):
            print("{config} not found run generate-config".format(config=args.config))
            exit(1)
        # Parse config
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).items()}
        working_folder = FolderHandler()
        working_folder.open_folder(os.path.join(config['output_dir'], "temp_trainModels"))

        samples = [SignalAlignSample(working_folder=working_folder, **s) for s in config["samples"]]
        # Train HMM transitions
        trainHMM(samples, working_folder, config, transitions=True)


if __name__ == "__main__":
    sys.exit(main())
