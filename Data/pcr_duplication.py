import sys
import argparse
from numpy import random
import math
from copy import deepcopy

def parse_args():
    defaults = {'number_of_cycles' : 6,
                'duplication_rate_per_cycle' : 0.75,
                'error_rate' : 1.0e-5,
                'random_seed' : 42,
                'pcr_product' : None}
    parser = argparse.ArgumentParser(
        description="Generates molecules")
    parser.add_argument("-i",
                        "--molecules",
                        type=str,
                        required=True,
                        help="Generated molecules in FASTA format")
    parser.add_argument("-c",
                        "--number-of-cycles",
                        type=int,
                        default=defaults['number_of_cycles'],
                        help="Number of PCR cycles (default: {})".format(defaults['number_of_cycles']))
    parser.add_argument("-d",
                        "--duplication-rate-per-cycle",
                        type=float,
                        default=defaults['duplication_rate_per_cycle'],
                        help="Proportion of molecules to be duplicated per cycle. (default: {})".format(defaults['duplication_rate_per_cycle']))
    parser.add_argument("-e",
                        "--error-rate",
                        type=float,
                        default=defaults['error_rate'],
                        help="PCR substituion error rate per duplicated base. (default: {})".format(defaults['error_rate']))
    parser.add_argument("-s",
                        "--random-seed",
                        type=int,
                        default=defaults['random_seed'],
                        help="NumPy random seed (default: {})".format(defaults['random_seed']))
    parser.add_argument("-o",
                        "--pcr-product",
                        type=str,
                        default=defaults['pcr_product'],
                        help="Output in FASTA format (default: stdout)")
    args = parser.parse_args()
    return args

def binary_search(length_index, position):
    start = 0
    end = len(length_index) - 1
    mid = (end - start) // 2 + start
    found = False
    while (not found):
        # print('OK:', start, mid, end, '\n',length_index, end='')
        if start >= end:
            print('WTF1')
            return
        if mid < start or end < mid:
            print('WTF3')
            return
        if length_index[mid] <= position and position < length_index[mid+1]:
            return mid, position - length_index[mid]
        elif position < length_index[mid]:
            end = mid
            mid = (end - start) // 2 + start
        elif length_index[mid + 1] <= position:
            start = mid + 1
            mid = (end - start) // 2 + start
        else:
            print('WTF2', end='')
            return

def pcr_cycle(molecules_list, duplication_rate, error_rate):
    duplication_size = int(len(molecules_list)*duplication_rate)
    duplication_sublist = [molecules_list[i] for i in random.choice(range(len(molecules_list)), size=duplication_size, replace=False)]
    duplicated_list = [(molecule_id, deepcopy(i)) for molecule_id,i in duplication_sublist]
    duplication_length = sum( (len(i[1]) for i in duplication_sublist) )
    length_index = [0]*(duplication_size+1)

    print(duplication_size, duplication_sublist, duplicated_list, duplication_length, length_index)
    for idx, molecule in enumerate(duplication_sublist):
        length_index[idx+1] = len(molecule[1]) + length_index[idx]
    duplication_length = length_index[-1]
    error_size = int(length_index[-1]*error_rate)
    print(length_index, error_size)
    for error_idx in random.choice(range(duplication_length), size=error_size, replace=False):
        idx, base = binary_search(length_index, error_idx)
        duplicated_list[idx][1][base] = chr(ord(duplicated_list[idx][1][base]) +1 )
    return duplicated_list

t = ['AAAA', 'CCCC', 'GGGG', 'TTTT']
t = [list(i) for i in t]

def main():
    args = parse_args()
    random.seed(args.random_seed)

    molecules_file = open(args.molecules)
    molecule_count = 0
    molecule_ids = []
    molecules = []
    line = molecules_file.readline()
    while(line != ''):
        molecule_ids.append(line.rstrip())
        line = molecules_file.readline()
        molecules.append(list(line.rstrip()))
    for i in molecules:
        print(i)

if __name__ == "__main__":
    main()
