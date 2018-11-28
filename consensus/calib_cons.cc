#include "spoa/spoa.hpp"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include <iostream>
#include <thread>
#include <functional>
#include <locale>

#define ALIGNMENT_TYPE_GLOBAL 1
#define SCORE_M 5
#define SCORE_X -3
#define SCORE_G -9
#define ASCII_SIZE 128
#define MSA_MAJORITY 0.5

int thread_count = 4;
int min_reads_per_cluster = 2;
std::string cluster_filename = "";
std::vector<std::string> fastq_filenames;
std::vector<std::string> output_filenames;
typedef uint32_t cluster_id_t;
typedef uint32_t read_id_t;


void print_help(){
    std::cout << "Calib Consensus: Generating consensus sequence from Calib clusters." << "\n";
    std::cout << "Usage: calib_cons [--PARAMETER VALUE]" << "\n";
    std::cout << "Example 1: calib_cons -t 8 -c input.cluster -q 1.fastq 2.fastq -o 1.out 2.out" << "\n";
    std::cout << "Example 2: calib_cons -q 1.fastq -q 2.fastq -o 1.out 2.out -c input.cluster" << "\n";
    std::cout << "Calib's paramters arguments:" << "\n";
    std::cout << "  -q  --fastq                    (type: space separated string list;\n";
    std::cout << "                                    REQUIRED paramter;\n";
    std::cout << "                                    can be set multiple times like in Example 2)\n";
    std::cout << "  -o  --output-prefix            (type: space separated string list;\n";
    std::cout << "                                    REQUIRED paramter;\n";
    std::cout << "                                    can be set multiple times like in Example 2;\n";
    std::cout << "                                    must be same size as fastq list)\n";
    std::cout << "  -c  --cluster                  (string;\n";
    std::cout << "                                    REQUIRED paramter)\n";
    std::cout << "  -t  --threads                  (positive integer;\n";
    std::cout << "                                    default: 4)\n";
    std::cout << "  -m  --min-reads-per-cluster    (positive integer;\n";
    std::cout << "                                    default: 2)\n";
    std::cout << "  -h  --help\n";
}

void parse_flags(int argc, char *argv[]){
    for (int i = 1; i < argc; i++) {
        std::string current_param(argv[i]);
        if (current_param == "-h" || current_param == "--help") {
            print_help();
            exit(0);
        }
        if (current_param == "-c" || current_param == "--cluster") {
            i++;
            cluster_filename = std::string(argv[i]);
            continue;
        }
        if (current_param == "-t" || current_param == "--threads") {
            i++;
            thread_count = atoi(argv[i]);
            continue;
        }
        if (current_param == "-m" || current_param == "--min-reads-per-cluster") {
            i++;
            min_reads_per_cluster = atoi(argv[i]);
            continue;
        }
        if (current_param == "-q" || current_param == "--fastq") {
            i++;
            for (; i < argc; i++) {
                std::string fastq_filename(argv[i]);
                if (fastq_filename[0] == '-') {
                    i--;
                    break;
                }
                fastq_filenames.push_back(fastq_filename);
            }
            continue;
        }
        if (current_param == "-o" || current_param == "--output-prefix") {
            i++;
            for (; i < argc; i++) {
                std::string output_filename(argv[i]);
                if (output_filename[0] == '-') {
                    i--;
                    break;
                }
                output_filenames.push_back(output_filename);
            }
            continue;
        }
        std::cout << "Unrecognized parameter, " << argv[i] << ", was passed.\n";
        print_help();
        exit(-1);
    }
    if (min_reads_per_cluster < 0 ) {
        std::cout << "Minimum reads per cluster ("<<min_reads_per_cluster<<") must be positive value." << '\n';
        print_help();
        exit(-1);
    }
    if (thread_count < 0 || thread_count > 16) {
        std::cout << "Thread count must be between 1 and 16." << '\n';
        print_help();
        exit(-1);
    }
    if (fastq_filenames.size() != output_filenames.size()) {
        std::cout << "Number of fastq files ("<< fastq_filenames.size() <<") must be equal number of output files ("<< output_filenames.size() <<")\n";
        print_help();
        exit(-1);
    }
    if (cluster_filename == "") {
        std::cout << "No cluster filename was passed.\n";
        print_help();
        exit(-1);
    }

}

void process_clusters(const std::vector<std::string>& read_to_sequence, const std::vector<std::vector<read_id_t> >& cluster_to_reads, std::string o_filename_prefix, size_t thread_id) {
    std::ofstream ofastq(o_filename_prefix + ".fastq" + std::to_string(thread_id));
    std::ofstream omsa(o_filename_prefix + ".msa" + std::to_string(thread_id));
    for (cluster_id_t cid = thread_id; cid < cluster_to_reads.size(); cid+=thread_count) {
        if (cluster_to_reads[cid].size() < min_reads_per_cluster) {
            continue;
        }
        auto alignment_engine = spoa::createAlignmentEngine(
            static_cast<spoa::AlignmentType>(ALIGNMENT_TYPE_GLOBAL), SCORE_M, SCORE_X, SCORE_G
        );
        auto graph = spoa::createGraph();
        std::stringstream header;
        header << cid << "\t";
        for (read_id_t rid : cluster_to_reads[cid]) {
            header << rid << ";";
            auto alignment = alignment_engine->align_sequence_with_graph(read_to_sequence[rid], graph);
            graph->add_alignment(alignment, read_to_sequence[rid]);
        }
        std::string consensus;
        std::string qual;
        std::vector<std::string> msa;
        graph->generate_multiple_sequence_alignment(msa);

        size_t profile_width = msa[0].size();
        size_t profile_height = msa.size();
        float profile[profile_width][ASCII_SIZE];
        for (int col = 0; col < profile_width; col++) {
            for (int row = 0; row < ASCII_SIZE; row++) {
                profile[col][row] = 0;
            }
        }
        for (const auto& it : msa) {
            for (int col = 0; col < profile_width; col++) {
                profile[col][it[col]]++;
            }
        }
        consensus.reserve(profile_width);
        for (int col = 0; col < profile_width; col++) {
            double majority_percentage = 0;
            bool is_gap = false;
            if        (profile[col]['A']/profile_height > MSA_MAJORITY) {
                consensus += 'A';
                majority_percentage = profile[col]['A']/profile_height;
            } else if (profile[col]['C']/profile_height > MSA_MAJORITY) {
                consensus += 'C';
                majority_percentage = profile[col]['C']/profile_height;
            } else if (profile[col]['G']/profile_height > MSA_MAJORITY) {
                consensus += 'G';
                majority_percentage = profile[col]['G']/profile_height;
            } else if (profile[col]['T']/profile_height > MSA_MAJORITY) {
                consensus += 'T';
                majority_percentage = profile[col]['T']/profile_height;
            } else if (profile[col]['-']/profile_height > MSA_MAJORITY) {
                is_gap = true;
            } else {
                consensus += 'N';
            }
            
            if (is_gap == false) {
                if (majority_percentage > 0.90) {
                    qual += 'K';
                } else if (majority_percentage > 0.70) {
                    qual += 'A';
                } else if (majority_percentage > 0.50) {
                    qual += '.';
                } else {
                    qual += '$';
                }
            }
        }
        ofastq << "@" << header.str() << '\n';
        ofastq << consensus << '\n';
        ofastq << '+' << '\n';
        ofastq << qual << '\n';

        omsa << "@" << header.str() << '\n';
        omsa << consensus << '\n';
        omsa << '+' << '\n';
        for (const auto& it: msa) {
            omsa << it << '\n';
        }
    }
}

void run_consensus(){
    std::vector<std::vector<read_id_t> > cluster_to_reads;
    size_t read_count = 0;
    std::string line_buffer;
    std::ifstream clusters_file;
    clusters_file.open(cluster_filename);
    // Get mapping of clusters to reads
    std::cout << "Reading cluster file: " << cluster_filename << '\n';
    while (getline(clusters_file, line_buffer)) {
        std::stringstream ss(line_buffer);
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> vstrings(begin, end);
        cluster_id_t cid = atoi(vstrings[0].c_str());
        read_id_t rid = atoi(vstrings[2].c_str());

        if (cluster_to_reads.size() <= cid) {
            cluster_to_reads.resize(cid+1);
        }
        cluster_to_reads[cid].push_back(rid);
        read_count++;
    }

    //Get read sequences from each FASTQ file, and pass it for MSA and output
    for (int i = 0; i < fastq_filenames.size(); i++) {
        std::cout << "Reading fastq file: " << fastq_filenames[i] << '\n';
        std::string ifastq_filename = fastq_filenames[i];
        std::cout << "Writing output files: " << output_filenames[i] << '\n';
        std::string o_filename_prefix = output_filenames[i];

        std::vector<std::string> read_to_sequence(read_count);
        std::ifstream ifastq;
        ifastq.open(ifastq_filename);
        read_id_t rid = 0;
        while (getline(ifastq, line_buffer)) {
            if (read_to_sequence.size() <= rid) {
                read_to_sequence.resize(rid+1);
            }
            getline(ifastq, read_to_sequence[rid]);
            getline(ifastq, line_buffer);
            getline(ifastq, line_buffer);
            rid++;
        }
        std::vector<std::thread> threads(thread_count);
        for (size_t thread_id = 0; thread_id < thread_count; thread_id++) {
            threads[thread_id] = std::thread(process_clusters, std::ref(read_to_sequence), std::ref(cluster_to_reads), o_filename_prefix, thread_id);
        }
        std::ofstream ofastq(o_filename_prefix + ".fastq");
        std::ofstream omsa(o_filename_prefix + ".msa");
        for (size_t thread_id = 0; thread_id < thread_count; thread_id++) {
            threads[thread_id].join();

            std::string fastq_t_filename = o_filename_prefix + ".fastq" + std::to_string(thread_id);
            std::ifstream ofastq_t(fastq_t_filename);
            ofastq << ofastq_t.rdbuf();
            remove(fastq_t_filename.c_str());

            std::string msa_t_filename = o_filename_prefix + ".msa" + std::to_string(thread_id);
            std::ifstream omsa_t(msa_t_filename);
            omsa << omsa_t.rdbuf();
            remove(msa_t_filename.c_str());
        }
    }
}

int main(int argc, char** argv) {
    parse_flags(argc, argv);
    run_consensus();
    return 0;
}
