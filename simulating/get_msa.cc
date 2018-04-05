#include <iostream>
#include <string>
#include <vector>

#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace std;
using namespace seqan;

int main () {
    string input_sequence;
    vector<string> sequences;
    while (getline (cin, input_sequence)) {
        sequences.push_back(input_sequence);
    }

    Align<Dna5String> align;
    resize(rows(align), sequences.size());
    for (int i = 0; i < sequences.size(); i++) {
        assignSource(row(align, i), sequences[i]);
    }

    globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));

    for (int i = 0; i < sequences.size(); i++) {
        cout << row(align,i) << "\n";
    }
    return 0;
}
