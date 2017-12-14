from pyfaidx import Fasta
import sys
import argparse
import random
import math


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generates molecules")
    parser.add_argument("-r",
                        "--reference",
                        type=str,
                        required=True,
                        help="Reference genome from which to generate molecules")
    parser.add_argument("-n",
                        "--number-of-molecules",
                        type=int,
                        default=1000,
                        help="Number of molecules. (default: 1000)")
    parser.add_argument("-u", "--molecule-size-mean",
                        type=int,
                        default=200,
                        help="Mean size of molecules (taken on standard dist. default: 200")
    parser.add_argument("-d",
                        "--molecule-size-standard-dev",
                        type=int,
                        default=25,
                        help="Molecule size stdev (default: 25)")
    parser.add_argument("-l",
                        "--min-molecule-size",
                        type=int,
                        default=150,
                        help="Any molecule size less than this will be adjusted to be equal to this by decreasing its start position. (default: 150)")
    parser.add_argument("-s",
                        "--random-seed",
                        type=int,
                        default=42,
                        help="NumPy random seed (default: 42)")
    parser.add_argument("-o", "--output-molecules", 
                        type=str,
                        default=None,
                        help="Output molecules/amplicons fasta file (default: stdout)")
    args = parser.parse_args()
    return args


def generate_molecule(genome_file_path,
                      output_file_path=sys.stdout,
                      molecule_size_mean=150,
                      molecule_size_standard_dev=20,
                      number_of_molecules=100,
                      random_seed=None,
                      min_molecule_size=150):

    if output_file_path:
        output_file = open(output_file_path, 'w+')
    else:
        output_file = sys.stdout
    random.seed(random_seed)
    genome = str(Fasta(genome_file_path)[0])
    genome_length = len(genome[0])
    for i in range(number_of_molecules):
        start = random.randrange(0, genome_length)
        end = min(
            [start + abs(math.floor(random.normalvariate(mu=molecule_size_mean, sigma=molecule_size_standard_dev))),
             genome_length])
        if end-start < min_molecule_size:
            start = end - min_molecule_size
        print(">{}:{}-{}".format(i, start, end), file=output_file)
        print(genome[start:end], file=output_file)


def main():
    args = parse_args()
    generate_molecule(genome_file_path=args.reference,
                      output_file_path=args.output_molecules,
                      molecule_size_mean=args.molecule_size_mean,
                      molecule_size_standard_dev=args.molecule_size_standard_dev,
                      number_of_molecules=args.number_of_molecules,
                      random_seed=args.random_seed,
                      min_molecule_size=args.min_molecule_size)


if __name__ == "__main__":
    main()
