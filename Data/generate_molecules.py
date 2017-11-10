from pyfaidx import Fasta
import sys
import argparse
import random
import math


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generates molecules/amplicons")
    parser.add_argument("-r", "--reference", type=str, help="Reference genome from which to generate molecules")
    parser.add_argument("-m", "--molecules", type=str, help="Output molecules/amplicons fasta file (default: stdout)",
                        default=None)
    parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: no seed)",
                        default=None)
    args = parser.parse_args()
    return args


def generate_molecule(genome_file_path,
                      output_file_path=sys.stdout,
                      molecule_size_mu=150,
                      molecule_size_sig=20,
                      number_of_molecules=100,
                      random_seed=None):

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
            [start + math.floor(random.normalvariate(mu=molecule_size_mu, sigma=molecule_size_sig)),
             genome_length])
        print('>'+str(i)+':'+str(start)+'-'+str(end), file=output_file)
        print(genome[0][start:end], file=output_file)


def main():
    args = parse_args()
    generate_molecule(args.reference, args.molecules, random_seed=args.random_seed)


if __name__ == "__main__":
    main()
