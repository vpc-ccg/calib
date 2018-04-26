#include <iostream>
#include <string>
#include <vector>

#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main (int argc, char *argv[]) {
    std::string input_sequence;
    if (argc != 3) {
        return -1;
    }
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    std::vector<std::string> sequences;
    while (getline (input, input_sequence)) {
        if (input_sequence.at(0) == '>') {
            sequences.clear();
            output << input_sequence << "\n";
        } else if (input_sequence.at(0) == '<') {
            Align<Dna5String> align;
            resize(rows(align), sequences.size());
            for (int i = 0; i < sequences.size(); i++) {
                assignSource(row(align, i), sequences[i]);
            }
            globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));
            for (int i = 0; i < sequences.size(); i++) {
                output << row(align,i) << "\n";
            }
            output << input_sequence << "\n";
        } else {
            sequences.push_back(input_sequence);
        }
    }

    return 0;
}
