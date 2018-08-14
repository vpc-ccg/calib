#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <seqan/align.h>
#include <seqan/graph_msa.h>

#define MSA_MAJORITY 0.5
#define ASCII_SIZE 128

std::string process_cluster(const std::vector<std::string> &sequences) {
    seqan::Align<seqan::Dna5String> align;
    seqan::resize(seqan::rows(align), sequences.size());
    for (int i = 0; i < sequences.size(); i++) {
        seqan::assignSource(seqan::row(align, i), sequences[i]);
    }
    std::string consensus = "";

    seqan::globalMsaAlignment(align, seqan::SimpleScore(5, -3, -9));

    std::vector<std::string> msa;
    for (int i = 0; i < sequences.size(); i++) {
        std::stringstream ss;
        ss << seqan::row(align,i);
        msa.push_back(ss.str());
    }
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
        if        (profile[col]['A']/profile_height > MSA_MAJORITY) {
            consensus += 'A';
        } else if (profile[col]['C']/profile_height > MSA_MAJORITY) {
            consensus += 'C';
        } else if (profile[col]['G']/profile_height > MSA_MAJORITY) {
            consensus += 'G';
        } else if (profile[col]['T']/profile_height > MSA_MAJORITY) {
            consensus += 'T';
        } else if (profile[col]['-']/profile_height > MSA_MAJORITY) {
            /* code */
        } else {
            consensus += 'N';
        }
    }
    return consensus;
}

int main (int argc, char *argv[]) {
    if (argc != 2) {
        return -1;
    }
    std::ifstream input(argv[1]);
    std::ofstream output(std::string(argv[1])+std::string(".seqan.fastq"));

    std::string input_sequence, consensus;
    std::vector<std::string> sequences;
    getline (input, consensus);
    while (getline (input, input_sequence)) {
        if (input_sequence.at(0) == '@') {
            consensus += '\t';
            consensus += process_cluster(sequences);
            output << consensus << "\n";
            sequences.clear();
            consensus = input_sequence;
        } else {
            sequences.push_back(input_sequence);
        }
    }
    consensus += '\t';
    consensus += process_cluster(sequences);
    output << consensus << "\n";
    sequences.clear();
    return 0;
}
