#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <sonLib.h>
#include <H5Fpublic.h>
#include <hdf5.h>

void usage() {
    fprintf(stderr, "\n\tsignalMachine - Align ONT ionic current to a reference sequence\n\n");
    fprintf(stderr, "--help: Display this super useful message and exit\n");
    fprintf(stderr, "-T: Template HMM model\n");
    fprintf(stderr, "-n: nucleotide fasta file location\n");
    fprintf(stderr, "-f: fast5 file location\n");
    fprintf(stderr, "-e: path to input event table in fast5\n");
    fprintf(stderr, "-E: path to output event table in fast5\n");
}



int main(int argc, char *argv[]) {
    char *templateModelFile = NULL;
    char *nucleotideFastaLocation= NULL;
    char *fast5Location = NULL;
    char *eventTableInput = NULL;
    char *eventTableOutput = NULL;

    int key;
    while (1) {
        static struct option long_options[] = {
                {"help",                    no_argument,        0,  'h'},
                {"templateModel",           required_argument,  0,  'T'},
                {"nucleotideFastaLocation", required_argument,  0,  'n'},
                {"fast5Location",           required_argument,  0,  'f'},
                {"eventTableInput",         required_argument,  0,  'e'},
                {"eventTableOutput",        required_argument,  0,  'E'},
                {0, 0, 0, 0} };

        int option_index = 0;

        key = getopt_long(argc, argv, "h:T:n:f:e:E:",
                          long_options, &option_index);

        if (key == -1) {
            //usage();
            break;
        }
        switch (key) {
            case 'h':
                usage();
                return 1;
            case 'T':
                templateModelFile = stString_copy(optarg);
                break;
            case 'n':
                nucleotideFastaLocation = stString_copy(optarg);
                break;
            case 'f':
                fast5Location = stString_copy(optarg);
                break;
            case 'e':
                eventTableInput = stString_copy(optarg);
                break;
            case 'E':
                eventTableOutput = stString_copy(optarg);
                break;
            default:
                usage();
                return 1;
        }
    }

    if (templateModelFile != NULL) printf("templateModelFile: %s\n", templateModelFile);
    if (nucleotideFastaLocation != NULL) printf("nucleotideFastaLocation: %s\n", nucleotideFastaLocation);
    if (fast5Location != NULL) printf("fast5Location: %s\n", fast5Location);
    if (eventTableInput != NULL) printf("eventTableInput: %s\n", eventTableInput);
    if (eventTableOutput != NULL) printf("eventTableOutput: %s\n", eventTableOutput);

    // open fast5 file
    hid_t file_id;
    if (!(file_id = H5Fopen(fast5Location, H5F_ACC_RDONLY , H5P_DEFAULT))) {
        st_errAbort("Could not open fast5 file: %s\n", fast5Location);
    }

    // sanity check file location
    if (!H5Lexists(file_id, eventTableInput, H5P_DEFAULT)) {
        st_errAbort("No entry in file at: %s\n", eventTableInput);
    }

    printf("\nfin.\n");
    return 0;
}
