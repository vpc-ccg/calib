#include <iostream>
#include <fstream>
//#include <string>
#include <stdio.h>
#include <map>

#include <cstdint> // include this header for uint64_t

using namespace std;
std::map<char, uint8_t> binarize;
std::map<uint8_t, char> charize ;

uint64_t reverse_complement(uint64_t kmer){
    uint64_t result = 0;
    //printf("%016lx\n", kmer);
    for (int i = 0; i < sizeof(uint64_t); i++){
        unsigned int which_byte = (sizeof(uint64_t) - i - 1);
        uint8_t  current_byte = kmer >> (8*which_byte) & 0xff;
        result |= ((uint64_t) current_byte ^ 0x03) << (8*i);
    }
    return result;
}

string bitvector_to_DNA(uint64_t kmer){
    char dna[8];
    for (int i = 0; i < sizeof(uint64_t); i++){
        dna[8-i-1] = (kmer >> (8*i)) & 0xff;
        dna[8-i-1] = charize[dna[8-i-1]];
    }
    return string(dna);
}


uint64_t minimizer(string seq, int start, int length){
    uint64_t current_kmer=0, current_rev_kmer=0, min_kmer= (uint64_t) - 1;

    for (unsigned int i = start; i < start + sizeof(uint64_t); i++){
        current_kmer <<= 8;
        current_rev_kmer >>= 8;
        uint8_t byte = binarize[seq.at(i)];
        current_kmer |= (uint64_t) byte;
        current_rev_kmer |= (uint64_t) (byte ^ 0x03 )<< 8*(sizeof(uint64_t) - 1);
    }

    for (int i = start + sizeof(uint64_t); i < start + length; i++){
        current_kmer <<= 8;
        current_rev_kmer >>= 8;
        uint8_t byte = binarize[seq.at(i)];
        current_kmer |= (uint64_t) byte;
        current_rev_kmer |= (uint64_t) (byte ^ 0x03 )<< 8*(sizeof(uint64_t) - 1);
        min_kmer = min_kmer < current_kmer ? min_kmer : current_kmer;
        min_kmer = min_kmer < current_rev_kmer ? min_kmer : current_rev_kmer;
    }
    return min_kmer;
}

int main(int argc, char** argv) {

    unsigned int barcode_length = 10;
    unsigned int k_size = 8;
    string out_filename = argv[2];
    out_filename = out_filename + ".out";

    ifstream fastq1;
    ifstream fastq2;
    fastq1.open (argv[1]);
    fastq2.open (argv[2]);
    ofstream output("myout");

    binarize['A'] = 0x00;
    binarize['C'] = 0x01;
    binarize['G'] = 0x02;
    binarize['T'] = 0x03;
    charize[0x00] = 'A';
    charize[0x01] = 'C';
    charize[0x02] = 'G';
    charize[0x03] = 'T';

    string r1, s1, q1, r2, s2, q2, trash;
    while (getline(fastq1, r1)) {
        getline(fastq1, s1);
        getline(fastq1, trash);
        getline(fastq1, q1);
        getline(fastq2, r2);
        getline(fastq2, s2);
        getline(fastq2, trash);
        getline(fastq2, q2);

        output << s1.substr(0, barcode_length) << "\t";
        output << s2.substr(0, barcode_length) << "\t";

        uint64_t min_kmer = minimizer(s1, barcode_length, 50);
        output << bitvector_to_DNA(min_kmer) << "\t";



        s1 = s1.substr(barcode_length);
        s2 = s2.substr(barcode_length);

    }

    return 0;
}


