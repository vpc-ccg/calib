import numpy as np
import argparse
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

def main():
    args = parse_args()
    np.random.seed(args.random_seed)
    molecules_file = open(args.molecules)
    molecule_names = list()
    molecule_seqs = list()
    molecule_cycles = list()
    molecule_parents = list()
    molecule_roots = list()
    line = molecules_file.readline()
    while len(line) > 0:
        molecule_names.append(line.rstrip())
        line = molecules_file.readline()
        molecule_seqs.append(line.rstrip())
        line = molecules_file.readline()
        molecule_cycles.append(0)
        molecule_parents.append(len(molecule_parents))
        molecule_roots.append(len(molecule_roots))
    mutations = dict(
        A=['C','G','T'],
        C=['A','G','T'],
        G=['C','A','T'],
        T=['C','G','A'],
        a=['C','G','T'],
        c=['A','G','T'],
        g=['C','A','T'],
        t=['C','G','A'],
        N=['A','C','G','T'],
        n=['A','C','G','T'],
    )
    for cycle in range(1, args.number_of_cycles+1):
        duplication_idxs_size = math.ceil(len(molecule_seqs)*args.duplication_rate_per_cycle)
        duplication_idxs_size = max(duplication_idxs_size, 2)

        duplication_idxs = list(np.random.choice(
                                len(molecule_seqs),
                                size=duplication_idxs_size,
                                replace=False))
        duplicated_molecules = list()
        duplicated_pos = list()
        for idx in duplication_idxs:
            duplicated_molecules.append(molecule_seqs[idx])
            for pos in range(len(duplicated_molecules[-1])):
                duplicated_pos.append((len(duplicated_molecules)-1, pos))
            molecule_parents.append(idx)
            molecule_roots.append(molecule_roots[idx])
            molecule_cycles.append(cycle)
        duplication_err_idxs_size = math.ceil(len(duplicated_pos)*args.error_rate)
        duplication_err_idxs_size = max(duplication_err_idxs_size, 2)

        duplication_err_idxs = list(np.random.choice(
                                len(duplicated_pos),
                                size=duplication_err_idxs_size,
                                replace=False))
        duplication_err_idxs = [duplicated_pos[idx] for idx in duplication_err_idxs]
        for (idx, pos) in duplication_err_idxs:
            base = duplicated_molecules[idx][pos]
            mutation = str(np.random.choice(mutations[base]))
            duplicated_molecules[idx] = '{}{}{}'.format(
                                            duplicated_molecules[idx][0:pos-1],
                                            mutation,
                                            duplicated_molecules[idx][pos+1:]
                                            )
        molecule_seqs.extend(duplicated_molecules)

    output_molecule_count = len(molecule_seqs)
    molecule_cycle_ancestry = [list()]*output_molecule_count
    for idx in range(output_molecule_count):
        ancestry = molecule_cycle_ancestry[molecule_parents[idx]].copy()
        ancestry.append(str(molecule_cycles[idx]))
        molecule_cycle_ancestry[idx] = ancestry
    output_molecules = list()
    for idx in range(output_molecule_count):
        output_molecules.append('{}:{}\n{}'.format(
            molecule_names[molecule_roots[idx]],
            '.'.join(molecule_cycle_ancestry[idx]),
            molecule_seqs[idx]
            )
        )
    np.random.shuffle(output_molecules)
    if args.pcr_product:
        output = open(args.pcr_product, 'w')
    else:
        output = sys.stdout
    for molecule in output_molecules:
        print(molecule, file=output)


if __name__ == "__main__":
    main()
