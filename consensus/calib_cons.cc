#include "spoa/spoa.hpp"
#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include <iostream>
#include <thread>
#include <mutex>
#include <functional>
#include <locale>

#define ALIGNMENT_TYPE_GLOBAL 1
#define SCORE_M 5
#define SCORE_X -3
#define SCORE_G -9
#define ASCII_SIZE 128
#define MSA_MAJORITY 0.5

size_t thread_count = 8;
std::mutex output_lock;

typedef uint32_t cluster_id_t;
typedef uint32_t read_id_t;

void process_clusters(const std::vector<std::string>& read_to_sequence, const std::vector<std::vector<read_id_t> >& cluster_to_reads, std::ofstream& ofastq, std::ofstream& omsa, size_t thread_id) {
    for (cluster_id_t cid = thread_id; cid < cluster_to_reads.size(); cid+=thread_count) {
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
        std::string spoa_consensus = graph->generate_consensus();
        std::stringstream calib_consensus;
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
        for (int col = 0; col < profile_width; col++) {
            if        (profile[col]['A']/profile_height > MSA_MAJORITY) {
                calib_consensus << 'A';
            } else if (profile[col]['C']/profile_height > MSA_MAJORITY) {
                calib_consensus << 'C';
            } else if (profile[col]['G']/profile_height > MSA_MAJORITY) {
                calib_consensus << 'G';
            } else if (profile[col]['T']/profile_height > MSA_MAJORITY) {
                calib_consensus << 'T';
            } else if (profile[col]['-']/profile_height > MSA_MAJORITY) {
                /* code */
            } else {
                calib_consensus << 'N';
            }
        }



        output_lock.lock();
        ofastq << ">" << header.str() << '\n';
        ofastq << calib_consensus.str() << '\n';
        ofastq << '+' << '\n';
        std::string qual(calib_consensus.str().size(), 'K');
        ofastq << qual << '\n';

        omsa << ">" << header.str() << '\n';
        omsa << calib_consensus.str() << '\n';
        omsa << spoa_consensus        << '\n';
        omsa << '+' << '\n';
        for (const auto& it: msa) {
            omsa << it << '\n';
        }
        output_lock.unlock();
    }
}

int main(int argc, char** argv) {
    std::vector<std::vector<read_id_t> > cluster_to_reads;
    size_t read_count = 0;

    std::string line_buffer;
    std::ifstream clusters_file;
    clusters_file.open(argv[1]);
    // Get mapping of clusters to reads
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
    for (int i = 2; i < argc; i++) {
        std::string filename = argv[i];

        std::vector<std::string> read_to_sequence(read_count);
        std::ifstream ifastq;
        ifastq.open(filename);
        read_id_t rid = 0;
        while (getline(ifastq, line_buffer)) {
            getline(ifastq, read_to_sequence[rid]);
            getline(ifastq, line_buffer);
            getline(ifastq, line_buffer);
            rid++;
        }
        std::ofstream ofastq(filename + ".consensus.fastq");
        std::ofstream omsa(filename + ".msa");

        std::vector<std::thread> threads(thread_count);
        for (size_t thread_id = 0; thread_id < thread_count; thread_id++) {
            threads[thread_id] = std::thread(process_clusters, std::ref(read_to_sequence), std::ref(cluster_to_reads), std::ref(ofastq), std::ref(omsa), thread_id);
        }
        for (size_t thread_id = 0; thread_id < thread_count; thread_id++) {
            threads[thread_id].join();
        }
    }
    return 0;
}
