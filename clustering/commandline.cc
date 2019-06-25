//
// Created by borabi on 19/12/17.
//

#include "commandline.h"

using namespace std;

#define READ_SIZE_SAMPLE_SIZE 10000

// Parameter definitions
string input_1 = "";
string input_2 = "";
string output_prefix = "";
bool silent = false;
bool no_sort = false;
int barcode_length_1 = -1;
int barcode_length_2 = -1;
int ignored_sequence_prefix_length = -1;
int minimizer_count = -1;
int error_tolerance = -1;
int minimizer_threshold = -1;
int kmer_size = -1;
int thread_count = -1;

void parse_flags(int argc, char *argv[]){
    int barcode_length = -1;
    for (int i = 1; i < argc; i++) {
        string current_param(argv[i]);
        if (current_param == "-h" || current_param == "--help") {
            print_help();
            exit(0);
        }
        if ((input_1 == "") && (current_param == "-f" || current_param == "--input-forward")) {
            input_1 = string(argv[i+1]);
            i++;
            continue;
        }
        if ((input_2 == "") && (current_param == "-r" || current_param == "--input-reverse")) {
            input_2 = string(argv[i+1]);
            i++;
            continue;
        }
        if ((output_prefix == "") && (current_param == "-o" || current_param == "--output-prefix")) {
            output_prefix = string(argv[i+1]);
            i++;
            continue;
        }
        if ((silent == false) && (current_param == "-s" || current_param == "--silent")) {
            silent = true;
            continue;
        }
        if ((no_sort == false) && (current_param == "-q" || current_param == "--no-sort")) {
            no_sort = true;
            continue;
        }
        if ((barcode_length == -1) && (barcode_length_1 == -1) && (barcode_length_2 == -1) && (current_param == "-l" || current_param == "--barcode-length")) {
            barcode_length_1 = atoi(argv[i+1]);
            barcode_length_2 = atoi(argv[i+1]);
            i++;
            continue;
        }
        if ((barcode_length_1 == -1) && (current_param == "-l1" || current_param == "--barcode-length-1")) {
            barcode_length_1 = atoi(argv[i+1]);
            i++;
            continue;
        }
        if ((barcode_length_2 == -1) && (current_param == "-l2" || current_param == "--barcode-length-2")) {
            barcode_length_2 = atoi(argv[i+1]);
            i++;
            continue;
        }
        if ((ignored_sequence_prefix_length == -1) && (current_param == "-p" || current_param == "--ignored-sequence-prefix-length")) {
            ignored_sequence_prefix_length = atoi(argv[i+1]);
            i++;
            continue;
        }
        if ((minimizer_count == -1) && (current_param == "-m" || current_param == "--minimizer-count")) {
            minimizer_count = atoi(argv[i+1]);
            i++;
            continue;
        }
        if ((kmer_size == -1) && (current_param == "-k" || current_param == "--kmer-size")) {
            kmer_size = atoi(argv[i+1]);
            i++;
            continue;
        }
        if ((error_tolerance == -1) && (current_param == "-e" || current_param == "--error-tolerance")) {
            error_tolerance = atoi(argv[i+1]);
            i++;
            continue;
        }
        if ((minimizer_threshold == -1) && (current_param == "-t" || current_param == "--minimizer-threshold")) {
            minimizer_threshold = atoi(argv[i+1]);
            i++;
            continue;
        }
        if ((thread_count == -1) && (current_param == "-c" || current_param == "--threads")) {
            thread_count = atoi(argv[i+1]);
            i++;
            continue;
        }

        cout << "Unrecognized parameter, repeated parameter, or parameter incompatible with previous parameters: " << current_param << "\n";
        print_help();
        exit(-1);
    }
    if ((barcode_length_1 + barcode_length_2 < 1) || (barcode_length_1 < 0 ) || (barcode_length_2 < 0)){
        cout << "Combined barcode lengths must be a positive integer and each mate barcode length must be non-negative! Note if both mates have the same barcode length you can use -l/--barcode-length parameter instead.\n";
        print_help();
        exit(-1);
    }
    if (input_1 == "" || input_2 == "" || output_prefix == "") {
        cout << "Missing input or output files parameters!\n";
        print_help();
        exit(-1);
    }
    if (ignored_sequence_prefix_length == -1) {
        ignored_sequence_prefix_length = 0;
    }
    if (thread_count == -1) {
        thread_count = 1;
    }
    if (error_tolerance == -1 && kmer_size == -1 &&  minimizer_count == -1 && minimizer_threshold == -1) {
        cout << "No error or minimizer parameters passed. Selecting parameters based on barcode and inferred read length\n";
        ifstream fastq1;
        fastq1.open (input_1);
        string line;
        size_t sample_read_count = 0;
        size_t total_reads_size = 0;
        while (getline(fastq1, line) && sample_read_count < READ_SIZE_SAMPLE_SIZE) {
            getline(fastq1, line);
            total_reads_size+=line.size();
            getline(fastq1, line);
            getline(fastq1, line);
            sample_read_count++;
        }
        size_t mean_read_size = total_reads_size/sample_read_count;
        int barcode_length = (barcode_length_1 + barcode_length_2)/2;
        if (barcode_length >= 1 && barcode_length <= 6) {
            if (mean_read_size >= 61 && mean_read_size <= 100) {
                error_tolerance     = 1;
                kmer_size           = 4;
                minimizer_count     = 6;
                minimizer_threshold = 2;
            }
            if (mean_read_size >= 101 && mean_read_size <= 150) {
                error_tolerance     = 1;
                kmer_size           = 8;
                minimizer_count     = 7;
                minimizer_threshold = 2;
            }
            if (mean_read_size >= 151 && mean_read_size <= 250) {
                error_tolerance     = 1;
                kmer_size           = 8;
                minimizer_count     = 7;
                minimizer_threshold = 2;
            }
        }
        if (barcode_length >= 7 && barcode_length <= 11) {
            if (mean_read_size >= 61 && mean_read_size <= 100) {
                error_tolerance     = 2;
                kmer_size           = 4;
                minimizer_count     = 7;
                minimizer_threshold = 3;
            }
            if (mean_read_size >= 101 && mean_read_size <= 150) {
                error_tolerance     = 2;
                kmer_size           = 8;
                minimizer_count     = 7;
                minimizer_threshold = 2;
            }
            if (mean_read_size >= 151 && mean_read_size <= 250) {
                error_tolerance     = 2;
                kmer_size           = 8;
                minimizer_count     = 7;
                minimizer_threshold = 2;
            }
        }
        if (barcode_length >= 12) {
            if (mean_read_size >= 61 && mean_read_size <= 100) {
                error_tolerance     = 2;
                kmer_size           = 4;
                minimizer_count     = 7;
                minimizer_threshold = 3;
            }
            if (mean_read_size >= 101 && mean_read_size <= 150) {
                error_tolerance     = 2;
                kmer_size           = 4;
                minimizer_count     = 7;
                minimizer_threshold = 3;
            }
            if (mean_read_size >= 151 && mean_read_size <= 250) {
                error_tolerance     = 2;
                kmer_size           = 8;
                minimizer_count     = 7;
                minimizer_threshold = 2;
            }
        }
        cout << "Inferred read length " << mean_read_size << " from sample of " << READ_SIZE_SAMPLE_SIZE << " reads\n";
        cout << "Selected paramters for (mean) barcode length " << barcode_length << " are:\n";
        cout << "\terror_tolerance\t" << error_tolerance << "\n";
        cout << "\tkmer_size\t" << kmer_size << "\n";
        cout << "\tminimizer_count\t" << minimizer_count << "\n";
        cout << "\tminimizer_threshold\t" << minimizer_threshold << "\n";
    }

    if (minimizer_count < 0 || error_tolerance < 0 || minimizer_threshold < 0 || kmer_size < 0) {
        cout << "Missing clustering error and minimizer parameters!\n";
        print_help();
        exit(-1);
    }
    if (thread_count < 1) {
        cout << "Number of threads must be between >= 1!\n";
        print_help();
        exit(-1);
    }
    if (minimizer_threshold > minimizer_count || minimizer_threshold < 0) {
        cout << "Minimizer threshold must be <= minimizer count\n";
        print_help();
        exit(-1);
    }
}

void print_flags(){
    cout << "Parameters:\n";
    cout << "\tinput_1:\t" << input_1 << "\n";
    cout << "\tinput_2:\t" << input_2 << "\n";
    cout << "\toutput_prefix:\t" << output_prefix << "\n";
    cout << "\tbarcode_length_1:\t" << barcode_length_1 << "\n";
    cout << "\tbarcode_length_2:\t" << barcode_length_2 << "\n";
    cout << "\tignored_sequence_prefix_length:\t" << ignored_sequence_prefix_length << "\n";
    cout << "\tminimizer_count:\t" << minimizer_count << "\n";
    cout << "\tkmer_size:\t" << kmer_size << "\n";
    cout << "\terror_tolerance:\t" << error_tolerance << "\n";
    cout << "\tminimizer_threshold:\t" << minimizer_threshold << "\n";
    cout << "\tthreads:\t" << thread_count << "\n";
    cout << "\n";

}

void print_help(){
    cout << "Calib: Clustering without alignment using LSH and MinHashing of barcoded reads" << "\n";
	cout << "Usage: calib [--PARAMETER VALUE]" << "\n";
	cout << "Example: calib -f R1.fastq -r R2.fastq -o my_out. -e 1 -l 8 -m 5 -t 2 -k 4 --silent" << "\n";
	cout << "Calib's paramters arguments:" << "\n";
    cout << "\t-f    --input-forward                 \t(type: string;   REQUIRED paramter)\n";
    cout << "\t-r    --input-reverse                 \t(type: string;   REQUIRED paramter)\n";
    cout << "\t-o    --output-prefix                 \t(type: string;   REQUIRED paramter)\n";
    cout << "\t-s    --silent                        \t(type: no value; default: unset)\n";
    cout << "\t-q    --no-sort                       \t(type: no value; default:  unset)\n";
    cout << "\t-l    --barcode-length                \t(type: int;      REQUIRED paramter unless -l1 and -l2 are provided)\n";
    cout << "\t-l1   --barcode-length-1              \t(type: int;      REQUIRED paramter unless -l is provided)\n";
    cout << "\t-l2   --barcode-length-2              \t(type: int;      REQUIRED paramter unless -l is provided)\n";
    cout << "\t-p    --ignored-sequence-prefix-length\t(type: int;      default: 0)\n";
    cout << "\t-m    --minimizer-count               \t(type: int;      default: Depends on observed read length;)\n";
    cout << "\t-k    --kmer-size                     \t(type: int;      default: Depends on observed read length;)\n";
    cout << "\t-e    --error-tolerance               \t(type: int;      default: Depends on observed read length;)\n";
    cout << "\t-t    --minimizer-threshold           \t(type: int;      default: Depends on observed read length;)\n";
    cout << "\t-c    --threads                       \t(type: int;      default: 1)\n";
    cout << "\t-h    --help\n";
}

// Memory use in MB
int get_memory_use(){
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            int i = strlen(line);
            const char* p = line;
            while (*p <'0' || *p > '9') p++;
            line[i-3] = '\0';
            result = atoi(p);
            break;
        }
    }
    fclose(file);
    return result/1024;
}
