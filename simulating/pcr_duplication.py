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
    root_molecule_seqs = list()
    molecule_cycles = list()
    molecule_parents = list()
    molecule_roots = list()
    molecule_mutations = list()
    line = molecules_file.readline()
    output_molecule_count = 0
    while len(line) > 0:
        molecule_names.append(line.rstrip())
        line = molecules_file.readline()
        root_molecule_seqs.append(list(line.rstrip()))
        line = molecules_file.readline()
        molecule_cycles.append(0)
        molecule_parents.append(len(molecule_parents))
        molecule_roots.append(len(molecule_roots))
        molecule_mutations.append(list())
        output_molecule_count += 1
    mutations = dict(
        # A=['C','G','T'],
        # C=['A','G','T'],
        # G=['C','A','T'],
        # T=['C','G','A'],
        # a=['C','G','T'],
        # c=['A','G','T'],
        # g=['C','A','T'],
        # t=['C','G','A'],
        # N=['A','C','G','T'],
        # n=['A','C','G','T'],
        A=['x'],
        C=['x'],
        G=['x'],
        T=['x'],
    )
    for cycle in range(1, args.number_of_cycles+1):
        print('Cycle {}'.format(cycle))
        parenting_count = math.ceil(output_molecule_count*args.duplication_rate_per_cycle)

        parenting_ids = list(np.random.choice(
                                output_molecule_count,
                                size=parenting_count,
                                replace=False))
        output_molecule_count += len(parenting_ids)

        cum_starts = list()
        cum_length = 0
        cum_starts.append(cum_length)
        for parent_id in parenting_ids:
            molecule_parents.append(parent_id)
            molecule_roots.append(molecule_roots[parent_id])
            molecule_cycles.append(cycle)
            molecule_mutations.append(molecule_mutations[parent_id].copy())
            cum_length += len(root_molecule_seqs[molecule_roots[-1]])
            cum_starts.append(cum_length)
        mutation_count = math.ceil(cum_length*args.error_rate)
        mutation_molecules = list(np.random.choice(
                                parenting_count,
                                size=mutation_count,
                                replace=True))
        for mol_id in mutation_molecules:
            mol_id = len(molecule_mutations) - parenting_count + mol_id
            mol_pos = np.random.randint(0, len(root_molecule_seqs[molecule_roots[mol_id]]))
            molecule_mutations[mol_id].append(mol_pos)

        # for cum_pos in mutation_cum_positions:
        #     guess_idx = int(cum_pos/cum_length*len(parenting_ids))
        #     search_iterator = None
        #     if cum_pos < cum_starts[guess_idx]:
        #         search_iterator = reversed(range(0, guess_idx))
        #     else:
        #         search_iterator = range(guess_idx, len(parenting_ids))
        #     for mol_id in search_iterator:
        #         if cum_pos >= cum_starts[mol_id] and cum_pos < cum_starts[mol_id+1]:
        #             mol_pos = cum_pos - cum_starts[mol_id]
        #             break
        #     mol_id = len(molecule_mutations) - parenting_count + mol_id
        #     molecule_mutations[mol_id].append(mol_pos)

    molecule_cycle_ancestry = [list()]*output_molecule_count
    for idx in range(output_molecule_count):
        ancestry = molecule_cycle_ancestry[molecule_parents[idx]].copy()
        ancestry.append(str(molecule_cycles[idx]))
        molecule_cycle_ancestry[idx] = ancestry
    output_molecules = list()
    for idx in range(output_molecule_count):
        molecule = root_molecule_seqs[molecule_roots[idx]].copy()
        for pos in molecule_mutations[idx]:
            molecule[pos] = str(np.random.choice(mutations[molecule[pos]]))
        output_molecules.append('{}:{}\t{}'.format(
                                molecule_names[molecule_roots[idx]],
                                '.'.join(molecule_cycle_ancestry[idx]),
                                ''.join(molecule),
                                ))
    np.random.shuffle(output_molecules)
    if args.pcr_product:
        output = open(args.pcr_product, 'w')
    else:
        output = sys.stdout
    for molecule in output_molecules:
        print(molecule, file=output)


if __name__ == "__main__":
    main()
