import numpy as np
import argparse
from numpy import random
import math
import sys

nucleotides = {'A', 'C', 'G', 'T'}


def parse_args():
    defaults = {'number_of_cycles': 6,
                'duplication_rate_per_cycle': 0.75,
                # Error rate for Taq is around 3e-5
                'error_rate': 3.0e-5,
                'random_seed': 42,
                'pcr_product': None}
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
                        help="PCR substitution error rate per duplicated base. (default: {})".format(defaults['error_rate']))
    parser.add_argument("-s",
                        "--random-seed",
                        type=int,
                        default=defaults['random_seed'],
                        help="Random seed (default: {})".format(defaults['random_seed']))
    parser.add_argument("-o",
                        "--pcr-product",
                        type=str,
                        default=defaults['pcr_product'],
                        help="Output in FASTA format (default: stdout)")
    args = parser.parse_args()
    return args


def flatten_dictionary(dictionary):
    flattened_dict = []
    for key, item in dictionary.items():
        for i in range(len(item)):
            new_tuple = (key, item[i])
            flattened_dict.append(new_tuple)
    return flattened_dict


def pcr_cycle(molecules, duplication_rate, error_rate, pcr_results, cycles_left=0):
    choices = flatten_dictionary(pcr_results)
    num_duplications = int(math.floor(len(choices) * duplication_rate))
    duplicate_choice_idxs = np.random.choice(len(choices), size=num_duplications, replace=False)
    print('PCR cycles left: ', cycles_left, file=sys.stderr)
    print('Select {} out of {} molecules'.format(num_duplications, len(choices)), file=sys.stderr)
    for i in duplicate_choice_idxs:
        molecule_id, existing_mutations = choices[i]
        molecule = molecules[molecule_id]
        error_choices = np.random.choice([True, False], p=[error_rate, 1-error_rate], size=len(molecule))
        error_idxs = np.where(error_choices)[0]
        mutations = []
        for error_idx in error_idxs:
            mutation_choices = nucleotides - set(molecule[error_idx])
            mutation = np.random.choice(list(mutation_choices), size=1)[0]
            mutations.append((error_idx, mutation))
        mutations = mutations + existing_mutations
        # print(mutations)
        pcr_results[molecule_id].append(mutations)
        # if len(molecules[molecule_id]) <= max(error_choices)
    if cycles_left == 1:
        return pcr_results
    else:
        return pcr_cycle(molecules, duplication_rate, error_rate, pcr_results, cycles_left-1)


def main():
    args = parse_args()
    random.seed(args.random_seed)
    np.random.seed(args.random_seed)
    molecules_file = open(args.molecules)
    pcr_results = {}
    molecules = {}
    line = 'init'
    while not line == '':
        line = molecules_file.readline()
        if line == '':
            break
        molecule_id = line.rstrip()
        molecule = molecules_file.readline().rstrip()
        molecules[molecule_id] = molecule
        # pcr_results is a list of lists of error indices
        pcr_results[molecule_id] = [[]]
    pcr_results = pcr_cycle(molecules,
                            args.duplication_rate_per_cycle,
                            args.error_rate,
                            pcr_results,
                            cycles_left=args.number_of_cycles)
    if args.pcr_product:
        output = open(args.pcr_product, 'w')
    else:
        output = sys.stdout
    final_molecule_count = 0
    for molecule_id, pcr_result in pcr_results.items():
        molecule_number = 0
        for mutations in pcr_result:
            molecule = molecules[molecule_id]
            for error_idx, mutation in mutations:
                molecule = molecule[0:error_idx] + mutation + molecule[error_idx+1:]
            print(">{}_{}_{}".format(final_molecule_count, molecule_id[1:], str(molecule_number)), file=output)
            molecule_number += 1
            final_molecule_count += 1
            print(molecule, file=output)


if __name__ == "__main__":
    main()
