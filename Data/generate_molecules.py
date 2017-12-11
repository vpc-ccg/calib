from pyfaidx import Fasta
import sys
import argparse
import random
import math


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generates molecules/amplicons")
    parser.add_argument("-r", "--reference", type=str, help="Reference genome from which to generate molecules",
                        required=True)
    parser.add_argument("-m", "--molecules", type=str, help="Output molecules/amplicons fasta file (default: stdout)",
                        default=None)
    parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: no seed)",
                        default=None)
    parser.add_argument("-u", "--molecule-size-mean", type=int,
                        help="Mean size of molecules (taken on standard dist. default: 150", default=150)
    parser.add_argument("-d", "--molecule-size-sig", type=int,
                        help="Molecule size stdev (default: 20)", default=20)
    parser.add_argument("-n", "--number-of-mol", type=int,
                        help="Number of molecules. (default: 100)", default=100)
    parser.add_argument("-l", "--min-molecule-size", type=int,
                        help="Any molecule size less than this will be adjusted to be equal to this by decreasing its start position. (default: 150)", default=150)
    args = parser.parse_args()
    return args


def generate_molecule(genome_file_path,
                      output_file_path=sys.stdout,
                      molecule_size_mu=150,
                      molecule_size_sig=20,
                      number_of_molecules=100,
                      random_seed=None,
                      min_molecule_size=150):

    if output_file_path:
        output_file = open(output_file_path, 'w+')
    else:
        output_file = sys.stdout
    random.seed(random_seed)
    genome = Fasta(genome_file_path)
    genome_length = len(genome[0])
    for i in range(number_of_molecules):
        start = random.randrange(0, genome_length)
        end = min(
            [start + abs(math.floor(random.normalvariate(mu=molecule_size_mu, sigma=molecule_size_sig))),
             genome_length])
        if end-start < min_molecule_size:
            start = end - min_molecule_size
        print(">{}:{}-{}".format(i, start, end), file=output_file)
        print(genome[0][start:end], file=output_file)


def main():
    args = parse_args()
    generate_molecule(genome_file_path=args.reference, 
                      output_file_path=args.molecules,
                      molecule_size_mu=args.molecule_size_mean,
                      molecule_size_sig=args.molecule_size_sig,
                      number_of_molecules=args.number_of_mol,
                      random_seed=args.random_seed,
                      min_molecule_size=args.min_molecule_size)


if __name__ == "__main__":
    main()
