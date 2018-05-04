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
    std::ofstream output(std::string(argv[2])+std::string(".stat"));
    std::ofstream output_msa(argv[2]);

    std::vector<std::string> sequences;
    while (getline (input, input_sequence)) {
        if (input_sequence.at(0) == '>') {
            sequences.clear();
            output << input_sequence.substr(1) << "\t";
            output_msa << input_sequence << "\n";
        } else if (input_sequence.at(0) == '<') {
            Align<Dna5String> align;
            resize(rows(align), sequences.size());
            for (int i = 0; i < sequences.size(); i++) {
                assignSource(row(align, i), sequences[i]);
            }
            globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));

            String<ProfileChar<Dna5> > profile;
            resize(profile, length(row(align, 0)));
            for (unsigned rowNo = 0; rowNo < sequences.size(); ++rowNo) {
                for (unsigned i = 0; i < length(row(align, rowNo)); ++i) {
                    profile[i].count[ordValue(getValue(row(align, rowNo), i))] += 1;
                }
            }

            int corrected_count = 0;
            int incorrigible_count = 0;
            for (unsigned i = 0; i < length(profile); ++i) {
                int max_idx = getMaxIndex(profile[i]);
                if (max_idx >= 0 && max_idx <= 4){
                    if (profile[i].count[max_idx] > sequences.size()/2) {
                        corrected_count += sequences.size() - profile[i].count[max_idx];
                    } else if (sequences.size() > 2) {
                        incorrigible_count += sequences.size();
                    }
                }
            }
            output << corrected_count << "\t" << incorrigible_count << "\n";
            for (int i = 0; i < sequences.size(); i++) {
                output_msa << row(align,i) << "\n";
            }
            output_msa << input_sequence << corrected_count << "\t" << incorrigible_count << "\n";

        } else {
            sequences.push_back(input_sequence);
        }
    }

    return 0;
}
