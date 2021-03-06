#!/usr/bin/env python3

from __future__ import print_function
import os
import sys
import csv
import numpy as np
import pysam
import subprocess
import h5py
import glob
from random import shuffle

import signalalign.utils.multithread as multithread
from signalalign import defaultModelFromVersion, parseFofn
from signalalign.nanoporeRead import NanoporeRead, NanoporeRead2D
from signalalign.utils.bwaWrapper import *
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.sequenceTools import fastaWrite, samtools_faidx_fasta, processReferenceFasta
from signalalign.mea_algorithm import mea_alignment_from_signal_align, match_events_with_signalalign


def create_signalAlignment_args(backward_reference=None, forward_reference=None, destination=None,
                                stateMachineType="threeState", in_templateHmm=None,
                                in_complementHmm=None, in_templateHdp=None, in_complementHdp=None, threshold=0.01,
                                diagonal_expansion=None,
                                constraint_trim=None, target_regions=None, degenerate=None, twoD_chemistry=False,
                                alignment_file=None, bwa_reference=None,
                                track_memory_usage=False, get_expectations=False, output_format='full', embed=False,
                                event_table=False,
                                check_for_temp_file_existance=True,
                                path_to_bin=''):
    """Create alignment arguments for SignalAlign. Parameters are explained in SignalAlignment"""
    alignment_args = {
        "backward_reference": backward_reference,
        "forward_reference": forward_reference,
        "destination": destination,
        "stateMachineType": stateMachineType,
        "in_templateHmm": in_templateHmm,
        "in_complementHmm": in_complementHmm,
        "in_templateHdp": in_templateHdp,
        "in_complementHdp": in_complementHdp,
        "threshold": threshold,
        "diagonal_expansion": diagonal_expansion,
        "constraint_trim": constraint_trim,
        "target_regions": target_regions,
        "degenerate": degenerate,
        "twoD_chemistry": twoD_chemistry,
        "alignment_file": alignment_file,
        "bwa_reference": bwa_reference,
        'track_memory_usage': track_memory_usage,
        'get_expectations': get_expectations,
        'output_format': output_format,
        'embed': embed,
        'event_table': event_table,
        'check_for_temp_file_existance': check_for_temp_file_existance,
        'path_to_bin': path_to_bin}

    return alignment_args


class SignalAlignment(object):
    def __init__(self,
                 in_fast5,
                 destination,
                 stateMachineType,
                 in_templateHmm,
                 in_complementHmm,
                 in_templateHdp,
                 in_complementHdp,
                 threshold,
                 diagonal_expansion,
                 constraint_trim,
                 degenerate,
                 forward_reference,
                 backward_reference=None,
                 # one of these needs to be set
                 alignment_file=None,
                 bwa_reference=None,
                 # reasonable defaults
                 twoD_chemistry=False,
                 target_regions=None,
                 output_format="full",
                 embed=False,
                 event_table=False,
                 check_for_temp_file_existance=True,
                 track_memory_usage=False,
                 get_expectations=False,
                 path_to_bin=''):
        self.in_fast5 = in_fast5  # fast5 file to align
        self.destination = destination  # place where the alignments go, should already exist
        self.stateMachineType = stateMachineType  # flag for signalMachine
        self.bwa_reference = bwa_reference  # path to reference sequence to generate guide alignment
        self.threshold = threshold  # min posterior probability to keep
        self.diagonal_expansion = diagonal_expansion  # alignment algorithm param
        self.constraint_trim = constraint_trim  # alignment algorithm param
        self.output_format = output_format  # smaller output files
        self.degenerate = degenerate  # set of nucleotides for degenerate characters
        self.twoD_chemistry = twoD_chemistry  # flag for 2D sequencing runs
        self.temp_folder = FolderHandler()  # object for holding temporary files (non-toil)
        self.read_name = self.in_fast5.split("/")[-1][:-6]  # get the name without the '.fast5'
        self.target_regions = target_regions
        self.output_formats = {"full": 0, "variantCaller": 1, "assignments": 2}
        self.embed = embed  # embed the output into the fast5 file
        self.event_table = event_table  # specify which event table to use to generate alignments
        self.backward_reference = backward_reference  # fasta path to backward reference if modified bases are used
        self.forward_reference = forward_reference  # fasta path to forward reference
        self.alignment_file = alignment_file  # guide aligments will be gotten from here if set
        self.check_for_temp_file_existance = check_for_temp_file_existance  # don't recreate if files exist
        self.track_memory_usage = track_memory_usage  # has the 'time' program append mem usage stats to output
        self.max_memory_usage_kb = None
        self.read_label = None
        self.get_expectations = get_expectations  # option to gather expectations of transitions and emissions
        self.path_to_signalMachine = os.path.join(path_to_bin, "signalMachine")  # path to signalMachine

        assert os.path.exists(self.path_to_signalMachine), "Path to signalMachine does not exist"
        assert self.bwa_reference is not None or self.alignment_file is not None, \
            "either 'bwa_reference' or 'alignment_file' argument is needed to generate cigar strings"

        if (in_templateHmm is not None) and os.path.isfile(in_templateHmm):
            self.in_templateHmm = in_templateHmm
        else:
            self.in_templateHmm = None
        if (in_complementHmm is not None) and os.path.isfile(in_complementHmm):
            self.in_complementHmm = in_complementHmm
        else:
            self.in_complementHmm = None

        # similarly for HDPs
        if (in_templateHdp is not None) and os.path.isfile(in_templateHdp):
            self.in_templateHdp = in_templateHdp
        else:
            self.in_templateHdp = None
        if (in_complementHdp is not None) and os.path.isfile(in_complementHdp):
            self.in_complementHdp = in_complementHdp
        else:
            self.in_complementHdp = None
        assert os.path.exists(self.destination), \
            "Destination path does not exist: {}".format(self.destination)

    def run(self):
        print("[SignalAlignment.run] INFO: Starting on {read}".format(read=self.in_fast5))
        if self.get_expectations:
            assert self.in_templateHmm is not None, "Need template HMM files for model training"
            if self.twoD_chemistry:
                assert self.in_complementHmm is not None, "Need compement HMM files for model training"
        if not os.path.isfile(self.in_fast5):
            print("[SignalAlignment.run] ERROR: Did not find .fast5 at{file}".format(file=self.in_fast5))
            return False

        # prep
        self.openTempFolder("tempFiles_%s" % self.read_name)
        if self.twoD_chemistry:
            npRead = NanoporeRead2D(fast_five_file=self.in_fast5, event_table=self.event_table, initialize=True)
        else:
            npRead = NanoporeRead(fast_five_file=self.in_fast5, event_table=self.event_table, initialize=True)
        #todo need to validate / generate events and nucleotide read

        # read label
        read_label = npRead.read_label  # use this to identify the read throughout
        self.read_label = read_label

        # nanopore read (event table, etc)
        npRead_ = self.addTempFilePath("temp_%s.npRead" % self.read_name)
        if not (self.check_for_temp_file_existance and os.path.isfile(npRead_)):
            # TODO is this totally fucked for RNA because of 3'-5' mapping?
            fH = open(npRead_, "w")
            ok = npRead.Write(out_file=fH, initialize=True)
            fH.close()
            if not ok:
                self.failStop("[SignalAlignment.run] File: %s did not pass initial checks" % self.read_name, npRead)
                return False

        # nucleotide read
        read_fasta_ = self.addTempFilePath("temp_seq_%s.fa" % read_label)
        ok = self.write_nucleotide_read(npRead, read_fasta_)
        if not ok:
            print("[SignalAlignment.run] Failed to write nucleotide read.  Continuing execution.")

        # alignment info
        cigar_file_ = self.addTempFilePath("temp_cigar_%s.txt" % read_label)
        temp_samfile_ = self.addTempFilePath("temp_sam_file_%s.sam" % read_label)
        strand = None
        reference_name = None
        if not (self.check_for_temp_file_existance and os.path.isfile(cigar_file_)):

            # need guide alignment to generate cigar file
            guide_alignment = None

            # get from alignment file
            if self.alignment_file is not None:
                guide_alignment = getGuideAlignmentFromAlignmentFile(self.alignment_file, read_name=read_label)
                if guide_alignment is None:
                    print("[SignalAlignment.run] read {} not found in {}".format(read_label, self.alignment_file))

            # get from bwa
            if guide_alignment is None and self.bwa_reference is not None:
                guide_alignment = generateGuideAlignment(reference_fasta=self.bwa_reference,
                                                         query=read_fasta_,
                                                         temp_sam_path=temp_samfile_,
                                                         target_regions=self.target_regions)
                if guide_alignment is None:
                    print("[SignalAlignment.run] read {} could not be aligned with BWA".format(read_label))

            # could not map
            if guide_alignment is None:
                self.failStop("[SignalAlignment.run] ERROR getting guide alignment", npRead)
                return False

            # ensure valid
            if not guide_alignment.validate():
                self.failStop("[SignalAlignment.run] ERROR invalid guide alignment", npRead)
                return False
            strand = guide_alignment.strand
            reference_name = guide_alignment.reference_name

            # write cigar to file
            cig_handle = open(cigar_file_, "w")
            cig_handle.write(guide_alignment.cigar + "\n")
            cig_handle.close()

        # otherwise, get strand from file
        else:
            strand, reference_name = getInfoFromCigarFile(cigar_file_)

        # add an indicator for the model being used
        if self.stateMachineType == "threeState":
            model_label = ".sm"
            stateMachineType_flag = ""
        elif self.stateMachineType == "threeStateHdp":
            model_label = ".sm3Hdp"
            stateMachineType_flag = "--sm3Hdp "
            if self.twoD_chemistry:
                assert (self.in_templateHdp is not None) and (self.in_complementHdp is not None), "Need to provide HDPs"
            else:
                assert self.in_templateHdp is not None, "Need to provide Template HDP"
        else:  # make invalid stateMachine control?
            model_label = ".sm"
            stateMachineType_flag = ""

        # next section makes the output file name with the format: /directory/for/files/file.model.orientation.tsv
        # forward strand
        if strand == "+":
            if self.output_format == "full":
                posteriors_file_path = os.path.join(self.destination, read_label + model_label + ".forward.tsv")
            elif self.output_format == "variantCaller":
                posteriors_file_path = os.path.join(self.destination, read_label + model_label + ".tsv")
            else:
                posteriors_file_path = os.path.join(self.destination, read_label + model_label + ".assignments.tsv")

        # backward strand
        elif strand == "-":
            if self.output_format == "full":
                posteriors_file_path = os.path.join(self.destination, read_label + model_label + ".backward.tsv")
            elif self.output_format == "variantCaller":
                posteriors_file_path = os.path.join(self.destination, read_label + model_label + ".tsv")
            else:
                posteriors_file_path = os.path.join(self.destination, read_label + model_label + ".assignments.tsv")

        # sanity check
        else:
            self.failStop("[SignalAlignment.run] ERROR Unexpected strand {}".format(strand), npRead)
            return False

        # flags

        # input (match) models
        if self.in_templateHmm is None:
            self.in_templateHmm = defaultModelFromVersion(strand="template", version=npRead.version)
        if self.twoD_chemistry and self.in_complementHmm is None:
            pop1_complement = npRead.complement_model_id == "complement_median68pA_pop1.model"
            self.in_complementHmm = defaultModelFromVersion(strand="complement", version=npRead.version,
                                                            pop1_complement=pop1_complement)

        assert self.in_templateHmm is not None
        if self.twoD_chemistry:
            if self.in_complementHmm is None:
                self.failStop("[SignalAlignment.run] ERROR Need to have complement HMM for 2D analysis", npRead)
                return False

        template_model_flag = "-T {} ".format(self.in_templateHmm)
        if self.twoD_chemistry:
            complement_model_flag = "-C {} ".format(self.in_complementHmm)
        else:
            complement_model_flag = ""

        print("[SignalAlignment.run] NOTICE: template model {t} complement model {c}"
              "".format(t=self.in_templateHmm, c=self.in_complementHmm))

        # reference sequences
        assert os.path.isfile(self.forward_reference)
        forward_ref_flag = "-f {f_ref} ".format(f_ref=self.forward_reference)
        if self.backward_reference:
            assert os.path.isfile(self.backward_reference)
            backward_ref_flag = "-b {b_ref} ".format(b_ref=self.backward_reference)
        else:
            backward_ref_flag = ""

        # input HDPs
        if (self.in_templateHdp is not None) or (self.in_complementHdp is not None):
            hdp_flags = "-v {tHdp_loc} ".format(tHdp_loc=self.in_templateHdp)
            if self.twoD_chemistry and self.in_complementHdp is not None:
                hdp_flags += "-w {cHdp_loc} ".format(cHdp_loc=self.in_complementHdp)
        else:
            hdp_flags = ""

        # threshold
        if self.threshold is not None:
            threshold_flag = "-D {threshold} ".format(threshold=self.threshold)
        else:
            threshold_flag = ""

        # diagonal expansion
        if self.diagonal_expansion is not None:
            diag_expansion_flag = "-x {expansion} ".format(expansion=self.diagonal_expansion)
        else:
            diag_expansion_flag = ""

        # constraint trim
        if self.constraint_trim is not None:
            trim_flag = "-m {trim} ".format(trim=self.constraint_trim)
        else:
            trim_flag = ""

        # output format
        if self.output_format not in list(self.output_formats.keys()):
            self.failStop("[SignalAlignment.run] ERROR illegal output format selected %s" % self.output_format)
            return False
        out_fmt = "-s {fmt} ".format(fmt=self.output_formats[self.output_format])

        # degenerate nucleotide information
        if self.degenerate is not None:
            degenerate_flag = "-o {} ".format(self.degenerate)
        else:
            degenerate_flag = ""

        # twoD flag
        if self.twoD_chemistry:
            twoD_flag = "--twoD"
        else:
            twoD_flag = ""

        # commands
        if self.get_expectations:
            template_expectations_file_path = os.path.join(self.destination, read_label + ".template.expectations.tsv")
            complement_expectations_file_path = os.path.join(self.destination,
                                                             read_label + ".complement.expectations.tsv")
            command = \
                "{vA} {td} {degen}{sparse}{model} -q {npRead} " \
                "{t_model}{c_model}{thresh}{expansion}{trim} {hdp}-L {readLabel} -p {cigarFile} " \
                "-t {templateExpectations} -c {complementExpectations} -n {seq_name} {f_ref_fa} {b_ref_fa}" \
                    .format(vA=self.path_to_signalMachine, model=stateMachineType_flag,
                            cigarFile=cigar_file_,
                            npRead=npRead_, readLabel=read_label, td=twoD_flag,
                            templateExpectations=template_expectations_file_path, hdp=hdp_flags,
                            complementExpectations=complement_expectations_file_path, t_model=template_model_flag,
                            c_model=complement_model_flag, thresh=threshold_flag, expansion=diag_expansion_flag,
                            trim=trim_flag, degen=degenerate_flag, sparse=out_fmt, seq_name=reference_name,
                            f_ref_fa=forward_ref_flag, b_ref_fa=backward_ref_flag)
        else:
            command = \
                "{vA} {td} {degen}{sparse}{model} -q {npRead} " \
                "{t_model}{c_model}{thresh}{expansion}{trim} -p {cigarFile} " \
                "-u {posteriors} {hdp}-L {readLabel} -n {seq_name} {f_ref_fa} {b_ref_fa}" \
                    .format(vA=self.path_to_signalMachine, model=stateMachineType_flag, sparse=out_fmt,
                            cigarFile=cigar_file_,
                            readLabel=read_label, npRead=npRead_, td=twoD_flag,
                            t_model=template_model_flag, c_model=complement_model_flag,
                            posteriors=posteriors_file_path, thresh=threshold_flag, expansion=diag_expansion_flag,
                            trim=trim_flag, hdp=hdp_flags, degen=degenerate_flag, seq_name=reference_name,
                            f_ref_fa=forward_ref_flag, b_ref_fa=backward_ref_flag)

        # run
        print("[SignalAlignment.run] running command: ", command, end="\n")
        try:
            command = command.split()

            if self.track_memory_usage:
                mem_command = ['/usr/bin/time', '-f', '\\nDEBUG_MAX_MEM:%M\\n']
                print("[SignalAlignment.run] Prepending command to track mem usage: {}".format(mem_command))
                mem_command.extend(command)
                command = mem_command

            output = subprocess.check_output(command, stderr=subprocess.STDOUT)
            output = str(output).split("\\n")
            for line in output:
                print("[SignalAlignment.run]    {}: {}".format(read_label, line))
                if line.startswith("DEBUG_MAX_MEM"):
                    self.max_memory_usage_kb = int(line.split(":")[1])

        except Exception as e:
            print("[SignalAlignment.run] exception ({}) running signalAlign: {}".format(type(e), e))
            raise e

        # save to fast5 file (if appropriate)
        if self.embed:
            print("[SignalAlignment.run] embedding into Fast5 ")

            data = self.read_in_signal_align_tsv(posteriors_file_path, file_type=self.output_format)
            npRead = NanoporeRead(fast_five_file=self.in_fast5, twoD=self.twoD_chemistry, event_table=self.event_table)
            npRead.Initialize(None)
            signal_align_path = npRead.get_latest_basecall_edition("/Analyses/SignalAlign_00{}", new=False)
            assert signal_align_path, "There is no path in Fast5 file: {}".format("/Analyses/SignalAlign_00{}")
            output_path = npRead._join_path(signal_align_path, self.output_format)
            npRead.write_data(data, output_path)

            # Todo add attributes to signalalign output
            if self.output_format == "full":
                print("[SignalAlignment.run] writing maximum expected alignment ")
                alignment = mea_alignment_from_signal_align(None, events=data)
                mae_path = npRead._join_path(signal_align_path, "MEA_alignment_labels")
                events = npRead.get_template_events()
                if events:
                    if strand == "-":
                        minus = True
                    else:
                        minus = False
                    labels = match_events_with_signalalign(sa_events=alignment,
                                                           event_detections=np.asanyarray(npRead.template_events),
                                                           minus=minus,
                                                           rna=npRead.is_read_rna())
                    npRead.write_data(labels, mae_path)
                    sam_string = str()
                    if os.path.isfile(temp_samfile_):
                        with open(temp_samfile_, 'r') as test:
                            for line in test:
                                sam_string += line
                    sam_path = npRead._join_path(signal_align_path, "sam")
                    # print(sam_string)
                    npRead.write_data(data=sam_string, location=sam_path, compression=None)

        # self.temp_folder.remove_folder()
        return True

    def write_nucleotide_read(self, nanopore_read, file_path):
        try:
            with open(file_path, "w") as read_file:
                # get appropriate read
                if self.twoD_chemistry:
                    # check for table to make 'assembled' 2D alignment table fasta with
                    if not nanopore_read.has2D_alignment_table:
                        nanopore_read.close()
                        return False
                    nucleotide_read = nanopore_read.alignment_table_sequence
                else:
                    nucleotide_read = nanopore_read.template_read

                # write read
                fastaWrite(fileHandleOrFile=read_file,
                           name=nanopore_read.read_label,
                           seq=nucleotide_read)

            return True
        except Exception as e:
            print('[SignalAlignment.write_nucleotide_read] {} exception: {}'.format(type(e), str(e)), file=sys.stderr)
            return False

    def openTempFolder(self, temp_dir):
        self.temp_folder.open_folder(os.path.join(self.destination, temp_dir))

    def addTempFilePath(self, path_to_add):
        return self.temp_folder.add_file_path(path_to_add)

    def failStop(self, message, nanopore_read=None):
        self.temp_folder.remove_folder()
        if nanopore_read is not None:
            nanopore_read.close()
        print(message)

    def read_in_signal_align_tsv(self, tsv_path, file_type):
        """Read in tsv file"""
        assert file_type in ("full", "assignments", "variantCaller")
        with open(tsv_path, 'r') as tsvin:
            if file_type == "full":
                dtype = [('contig', 'S10'), ('reference_index', int),
                         ('reference_kmer', 'S5'), ('read_file', 'S57'),
                         ('strand', 'S1'), ('event_index', int),
                         ('event_mean', float), ('event_noise', float),
                         ('event_duration', float), ('aligned_kmer', 'S5'),
                         ('scaled_mean_current', float), ('scaled_noise', float),
                         ('posterior_probability', float), ('descaled_event_mean', float),
                         ('ont_model_mean', float), ('path_kmer', 'S5')]
            elif file_type == "assignments":
                dtype = [('k-mer', 'S10'), ('read_file', 'S57'),
                         ('descaled_event_mean', float), ('posterior_probability', float)]

            else:
                dtype = [('event_index', int), ('reference_position', int),
                         ('base', 'S6'), ('posterior_probability', float), ('strand', 'S1'),
                         ('forward_mapped', int), ('read_file', 'S57')]

            event_table = np.loadtxt(tsvin, dtype=dtype)

            def remove_field_name(a, name):
                names = list(a.dtype.names)
                if name in names:
                    names.remove(name)
                b = a[names]
                return b

            event_table = remove_field_name(event_table, "read_file")

        return event_table


def signal_alignment_service(work_queue, done_queue, service_name="signal_alignment"):
    # prep
    total_handled = 0
    failure_count = 0
    mem_usages = list()

    # catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                alignment = SignalAlignment(**f)
                success = alignment.run()
                if alignment.max_memory_usage_kb is not None:
                    mem_usages.append(alignment.max_memory_usage_kb)
                if not success: failure_count += 1
            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, multithread.current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, multithread.current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, multithread.current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(multithread.TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(multithread.FAILURE_KEY, failure_count))
        if len(mem_usages) > 0:
            done_queue.put("{}:{}".format(multithread.MEM_USAGE_KEY, ",".join(map(str, mem_usages))))


def multithread_signal_alignment(signal_align_arguments, fast5_locations, worker_count, forward_reference=None):
    """Multiprocess SignalAlignment for a list of fast5 files given a set of alignment arguments.

    :param signal_align_arguments: signalAlignment arguments besides 'in_fast5'
    :param fast5_locations: paths to fast5 files
    :param worker_count: number of workers
    :param forward_reference: path to forward reference for signalAlign alignment
    """
    # don't modify the signal_align_arguments
    signal_align_arguments = dict(**signal_align_arguments)
    if not forward_reference:
        assert "forward_reference" in signal_align_arguments, "Must specify forward_reference path"
        forward_reference = signal_align_arguments['forward_reference']

    # Samtools Index reference files for quick access
    samtools_faidx_fasta(forward_reference, log="multithread_signal_alignment")
    if "backward_reference" in signal_align_arguments and signal_align_arguments["backward_reference"]:
        samtools_faidx_fasta(signal_align_arguments["backward_reference"], log="multithread_signal_alignment")

    # bwa_reference and alignment_file must be specified
    assert ('bwa_reference'  in signal_align_arguments and signal_align_arguments["bwa_reference"]  is not None) or \
           ('alignment_file' in signal_align_arguments and signal_align_arguments['alignment_file'] is not None), \
        "Must specify bwa_reference or alignment_file"

    # if we didn't get an alignment file then check for index file
    if 'alignment_file' not in signal_align_arguments or not signal_align_arguments['alignment_file']:
        # ensure alignments can be generated (either from bwa on the reference or by an alignment file)
        bwa_reference = buildBwaIndex(signal_align_arguments["bwa_reference"],
                                      os.path.dirname(signal_align_arguments["bwa_reference"]),
                                      log='multithread_signal_alignment')
        signal_align_arguments["bwa_reference"] = bwa_reference
        assert os.path.exists(bwa_reference + ".bwt"), "Error creating BWA index for: {}".format(bwa_reference)

    # ensure required arguments are in signal_align_argments
    required_arguments = {'destination', 'stateMachineType', 'in_templateHmm', 'in_complementHmm',
                          'in_templateHdp', 'in_complementHdp', 'threshold', 'diagonal_expansion', 'constraint_trim',
                          'degenerate', 'forward_reference'}
    optional_arguments = {'backward_reference', 'alignment_file', 'bwa_reference', 'twoD_chemistry',
                          'target_regions', 'output_format', 'embed', 'event_table', 'check_for_temp_file_existance',
                          'track_memory_usage', 'get_expectations', 'path_to_bin'}
    missing_arguments = list(filter(lambda x: x not in signal_align_arguments.keys(), required_arguments))
    unexpected_arguments = list(filter(lambda x: x not in required_arguments and x not in optional_arguments,
                                       signal_align_arguments.keys()))
    assert len(missing_arguments) == 0 and len(unexpected_arguments) == 0, \
        "Invalid arguments to signal_align.  Missing: {}, Invalid: {}".format(missing_arguments, unexpected_arguments)
    assert os.path.exists(signal_align_arguments['destination']), \
        "Destination path does not exist: {}".format(signal_align_arguments['destination'])
    # run the signal_align_service
    print("[multithread_signal_alignment] running signal_alignment on {} fast5s with {} workers".format(
        len(fast5_locations), worker_count))
    total, failure, messages = multithread.run_service(
        signal_alignment_service, fast5_locations, signal_align_arguments, 'in_fast5', worker_count)

    # report memory usage
    memory_stats = list()
    for message in messages:
        if message.startswith(multithread.MEM_USAGE_KEY):
            memory_stats.extend(map(int, message.split(":")[1].split(",")))
    if len(memory_stats) > 0:
        kb_to_gb = lambda x: float(x) / (1 << 20)
        print("[multithread_signal_alignment] memory avg: %3f Gb" % (kb_to_gb(np.mean(memory_stats))))
        print("[multithread_signal_alignment] memory std: %3f Gb" % (kb_to_gb(np.std(memory_stats))))
        print("[multithread_signal_alignment] memory max: %3f Gb" % (kb_to_gb(max(memory_stats))))

    # fin
    print("[multithread_signal_alignment] fin")
    return [x for x in glob.glob(os.path.join(signal_align_arguments["destination"], "*.tsv")) if
            os.stat(x).st_size != 0]


def create_sa_sample_args(fofns=[], fast5_dirs=[], positions_file=None, motifs=None, alignment_file=None,
                          bwa_reference=None, fw_reference=None, bw_reference=None, name=None,
                          number_of_kmer_assignments=10, probability_threshold=0.8, kmers_from_reference=False):
    """Create sample arguments for SignalAlignSample. Parameters are explained in SignalAlignmentSample"""
    sample_args = {
        "fofns": fofns,
        "fast5_dirs": fast5_dirs,
        "positions_file": positions_file,
        "motifs": motifs,
        "bwa_reference": bwa_reference,
        "fw_reference": fw_reference,
        "bw_reference": bw_reference,
        "name": name,
        "number_of_kmer_assignments": number_of_kmer_assignments,
        "probability_threshold": probability_threshold,
        "kmers_from_reference": kmers_from_reference,
        'alignment_file': alignment_file
    }
    return sample_args


class SignalAlignSample(object):
    def __init__(self, working_folder, fofns, fast5_dirs, positions_file, motifs, bwa_reference, fw_reference,
                 bw_reference, name, number_of_kmer_assignments, probability_threshold, kmers_from_reference,
                 alignment_file):
        """Prepare sample for processing via signalAlign.

        :param working_folder: FolderHandler() object with a working directory already created
        :param fofns: list of fofns aka 'File Of File Names'
        :param fast5_dirs: list of fast5 directories
        :param positions_file: a path to tsv file. format: [contig  position    strand  change_from change_to]
        :param motifs: list of motifs to process canonical guide reference fasta
        :param bwa_reference: canoncial reference genome fasta. (preferably bwa indexed already)
        :param fw_reference: path to the forward signalAlign reference (edited if looking for ambiguous positions or motifs)
        :param bw_reference: path to the backward signalAlign reference (edited if looking for ambiguous positions or motifs)
        :param name: name of sample
        :param alignment_file: path to bam file to get alignments

        ###### Sample specific HDP Training Parameters #######
        :param number_of_kmer_assignments: max number of assignments for each kmer
        :param probability_threshold: minimum probability required to use for training
        :param kmers_from_reference: extract training kmers from reference sequence.

        """
        self.motifs = motifs
        self.bwa_reference = bwa_reference
        self.alignment_file = alignment_file
        self.fw_fasta_path = fw_reference
        self.bw_fasta_path = bw_reference
        self.positions_file = positions_file
        self.fast5_dirs = fast5_dirs
        self.fofns = fofns
        self.name = name
        self.number_of_kmer_assignments = number_of_kmer_assignments
        self.probability_threshold = probability_threshold
        self.kmers_from_reference = kmers_from_reference
        self.working_folder = working_folder

        assert self.name is not None, "Must specify a name for your sample. name: {}".format(self.name)
        assert isinstance(self.fast5_dirs, list), "fast5_dirs needs to be a list. fast5_dirs: {}".format(
            self.fast5_dirs)
        assert isinstance(self.fofns, list), "fofns needs to be a list. fofns: {}".format(self.fofns)
        assert len(self.fast5_dirs) + len(self.fofns) >= 1, "must specifcy at least one fofn and one fast5_dir. " \
                                                            "fofns: {}, fast5_dirs: {}".format(self.fofns,
                                                                                               self.fast5_dirs)
        # container for analysis files
        self.analysis_files = None
        self.files = []
        self._find_fast5_files()
        self.process_references()

    def _find_fast5_files(self):
        """Get all fasta paths via fofn.txt files and fast5_dirs"""
        # gather all files in fast5 directories
        for fast5_dir in self.fast5_dirs:
            assert os.path.isdir(fast5_dir), "Fast5 directory does not exist. dir: {}".format(fast5_dir)
            self.files.extend(
                [os.path.abspath(os.path.join(fast5_dir, x)) for x in os.listdir(fast5_dir) if x.endswith(".fast5")])
        # gather files from fofn.txt files
        for fofn in self.fofns:
            self.files.extend(parseFofn(fofn))
        assert len(self.files) > 0, "Found no files in fast5_dirs or fofns. files: " \
                                    "{}, fofns: {}, fast5_dirs: {}".format(self.files, self.fofns, self.fast5_dirs)

    def getFiles(self):
        return self.files

    def getReferences(self):
        return self.fw_fasta_path, self.bw_fasta_path

    def process_references(self):
        """Process a set of Fast5 files. Creates edited reference sequences if needed"""
        if self.fw_fasta_path is None:
            assert os.path.isfile(
                self.bwa_reference), "Must specify a bwa_reference in order to create signalAlignments. {}" \
                                     "".format(self.bwa_reference)
            self.fw_fasta_path, self.bw_fasta_path = processReferenceFasta(fasta=self.bwa_reference,
                                                                           work_folder=self.working_folder,
                                                                           motifs=self.motifs,
                                                                           positions_file=self.positions_file)


# TODO use Fast5 object
def get_2d_length(fast5, verbose=False):
    """Get 2d Read length. Searches only in one location

   :param fast5: path to fast5 file
   :param verbose: bool option to print update info
    """
    read = h5py.File(fast5, 'r')
    twoD_read_sequence_address = "/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
    if not (twoD_read_sequence_address in read):
        if verbose:
            print("This read didn't have a 2D read", fast5, end='\n', file=sys.stderr)
        read.close()
        return 0
    else:
        read_length = len(read[twoD_read_sequence_address][()].split()[2])
        if verbose:
            print("read %s has %s bases" % (fast5, read_length))
        read.close()
        return read_length


# TODO use Fast5 object
def get_1d_length(fast5, verbose=False):
    """Get 1D Read length. Searches only in one location

    :param fast5: path to fast5 file
    :param verbose: bool option to print update info
    """
    read = h5py.File(fast5, "r")
    template_fastq_address = "/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    if not (template_fastq_address in read):
        if verbose:
            print("Read %s has not been basecalled" % fast5)
        read.close()
        return 0
    else:
        read_length = len(read[template_fastq_address][()].split()[2])
        if verbose:
            print("read %s has %s bases" % (fast5, read_length))
        read.close()
        return read_length


def trim_num_files_in_sample(sample, max_bases, twoD, verbose=True):
    """Trim the number of fast5 files based on sequence length
    :param sample: "AbstractSamples" samples
    :param max_bases: max number of nucelotides to analyze
    :param twoD: boolean option to include twoD length
    :param verbose: boolean option to report trimming stats
    """
    shuffle(sample.getFiles())
    total_amount = 0
    file_count = 0
    get_seq_len_fcn = get_2d_length if twoD else get_1d_length
    # loop over files and add them to training list, break when we have enough bases to complete a batch
    # collect paths to fast5 files
    list_of_fast5s = []
    for f in sample.getFiles():
        total_amount += get_seq_len_fcn(f)
        if total_amount >= max_bases:
            if len(list_of_fast5s) == 0:
                total_amount = 0
                continue
            else:
                break
        # yield f, sample.fw_fasta_path, sample.bw_fasta_path
        list_of_fast5s.append(f)
        # training_files.append((f, sample.fw_fasta_path, sample.bw_fasta_path))
        file_count += 1
        if verbose:
            print("[Trim_Sample] Culled {file_count} training files, for {bases} from {sample}."
                  .format(file_count=file_count, bases=total_amount, sample=sample.name, end="\n", file=sys.stderr))
    assert len(list_of_fast5s) > 0, "All fast5 files have sequence lengths greater than max_bases: {}".format(max_bases)
    return list_of_fast5s


def multithread_signal_alignment_samples(samples, signal_align_arguments, worker_count, trim=None):
    """Multiprocess SignalAlignment for a list of fast5 files given a set of alignment arguments.

    :param samples: list of "process_sample" samples
    :param signal_align_arguments: signalAlignment arguments besides 'in_fast5'
    :param worker_count: number of workers
    :param trim: number of bases to analyze for each sample
    """
    original_destination = signal_align_arguments["destination"]
    names = [sample.name for sample in samples]
    assert len(names) == len(set(names)), "Cannot have same name for multiple samples: sample.name {}" \
                                          "".format(names)
    for sample in samples:
        # process sample
        if trim:
            assert type(trim) is int, "Trim must be an integer"
            list_of_fast5s = trim_num_files_in_sample(sample, trim, signal_align_arguments["twoD_chemistry"],
                                                      verbose=True)
        else:
            list_of_fast5s = sample.getFiles()
        # correct signal align arguments
        signal_align_arguments["alignment_file"] = sample.alignment_file
        signal_align_arguments["bwa_reference"] = sample.bwa_reference
        signal_align_arguments["backward_reference"] = sample.bw_fasta_path
        signal_align_arguments["forward_reference"] = sample.fw_fasta_path
        signal_align_arguments["destination"] = os.path.join(original_destination, sample.name)
        if not os.path.exists(signal_align_arguments["destination"]):
            os.mkdir(signal_align_arguments["destination"])
        # run signal align
        output_files = multithread_signal_alignment(signal_align_arguments, list_of_fast5s, worker_count)
        sample.analysis_files = output_files

    return samples
