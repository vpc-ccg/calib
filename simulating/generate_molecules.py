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
    parser.add_argument("-b",
                        "--bed",
                        type=str,
                        help="bed file that determines targeted regions of genome")
    args = parser.parse_args()
    return args

def generate_molecules_from_bed(genome_file_path,
                                bed_file_path,
                                molecule_size_mean,
                                molecule_size_standard_dev,
                                number_of_molecules,
                                random_seed,
                                min_molecule_size,
                                output_file_path=None):
    random.seed(random_seed)
    genome_obj = Fasta(genome_file_path)
    genome = dict()
    for key in genome_obj.keys():
        genome[key] = str(genome_obj[key])
    del genome_obj
    bed_positions = set()
    for line in open(bed_file_path).readlines()[1:]:
        line = line.rstrip().split('\t')
        chr = line[0]
        start = int(line[1])
        end = int(line[2])
        gene = line[3]
        for idx in range(start, end + 1):
            bed_positions.add((chr, idx, gene))
    bed_positions = list(bed_positions)

    molecules = list()
    while len(molecules) < number_of_molecules:
        mol_len = math.floor(random.normalvariate(mu=molecule_size_mean, sigma=molecule_size_standard_dev))
        if mol_len < min_molecule_size:
            continue
        chr, pos, gene = random.choice(bed_positions)
        start = pos - math.floor(random.normalvariate(mu=mol_len/2, sigma=15))
        end   = pos + mol_len
        if start <= 0 or end > len(genome[chr]):
            continue
        molecules.append((chr, gene, start, end))

    if output_file_path:
        output_file = open(output_file_path, 'w+')
    else:
        output_file = sys.stdout

    molecules.sort()
    for idx, molecule in enumerate(molecules):
        chr   = molecule[0]
        gene  = molecule[1]
        start = molecule[2]
        end   = molecule[3]
        seq = str(genome[chr][start: end]).upper()

        print(">{}:{}:{}:{}:{}".format(idx, chr, gene, start, end), file=output_file)
        print(seq, file=output_file)

def main():
    args = parse_args()
    if args.bed:
        generate_molecules_from_bed(genome_file_path=args.reference,
                                    bed_file_path=args.bed,
                                    molecule_size_mean=args.molecule_size_mean,
                                    molecule_size_standard_dev=args.molecule_size_standard_dev,
                                    number_of_molecules=args.number_of_molecules,
                                    random_seed=args.random_seed,
                                    min_molecule_size=args.min_molecule_size,
                                    output_file_path=args.output_molecules)


if __name__ == "__main__":
    main()
