#include "spoa/spoa.hpp"
#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include <unordered_map>
#include <iostream>

#define ALIGNMENT_TYPE_GLOBAL 1
#define SCORE_M 5
#define SCORE_X -4
#define SCORE_G -8

typedef uint32_t cluster_id_t;
typedef uint32_t read_id_t;
int main(int argc, char** argv) {
    std::vector<std::vector<read_id_t> > cluster_to_reads;
    std::vector<cluster_id_t> read_to_cluster;

    std::string line_buffer;
    std::ifstream clusters_file;
    clusters_file.open(argv[1]);
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

        if (read_to_cluster.size() <= rid) {
            read_to_cluster.resize(rid+1);
        }
        read_to_cluster[rid] = cid;
    }

    for (int i = 2; i < argc; i++) {
        std::string filename = argv[i];

        std::vector<std::string> read_to_sequence(read_to_cluster.size());
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

        for (cluster_id_t cid = 0; cid < cluster_to_reads.size(); cid++) {
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
            std::string consensus = graph->generate_consensus();
            std::vector<std::string> msa;
            graph->generate_multiple_sequence_alignment(msa);

            ofastq << ">" << header.str() << '\n';
            ofastq << consensus << '\n';
            ofastq << '+' << '\n';
            std::string qual(consensus.size(), 'K');
            ofastq << qual << '\n';

            omsa << ">" << header.str() << '\n';
            omsa << consensus << '\n';
            omsa << '+' << '\n';
            for (const auto& it: msa) {
                omsa << it << '\n';
            }
        }
    }
    return 0;
}
