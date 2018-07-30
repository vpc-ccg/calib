from pyfaidx import Fasta
import pandas as pd
import numpy as np
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


def generate_random_molecule(mu, sigma, genome_length, chromosome_ends, min_molecule_size):
    start = random.randrange(0, genome_length)
    chromosome_indices = np.array(list(chromosome_ends.values()))
    distances = chromosome_indices - start
    chrom_end = np.min(distances[distances > 0])
    chromosome = np.where(distances == chrom_end)[0][0]

    end = min(
        [start + abs(math.floor(random.normalvariate(mu=mu, sigma=sigma))),
         chrom_end])
    while end - start < 0:
        end = min(
            [start + abs(math.floor(random.normalvariate(mu=mu, sigma=sigma))),
             genome_length])
    if end - start < min_molecule_size:
        start = end - min_molecule_size
    return chromosome, start, end


def process_genome(genome):
    # genome is a pyfaidx.Fasta object
    chromosome_strings = [str(genome[chrom_id]) for chrom_id in genome.keys()]
    chromosome_ids = list(genome.keys())
    hg38_chromosomes = [str(i) for i in range(1, 23)] + ['Y', 'X']
    genome_str = "".join(chromosome_strings)
    chromosome_indices = np.array([len(chrom) for chrom in chromosome_strings])
    running_sum = 0
    chromosome_ends = {}
    for chrom_id, chrom_length in zip(chromosome_ids, chromosome_indices):
        if chrom_id in hg38_chromosomes:
            chrom_id = 'chr' + chrom_id
        running_sum += chrom_length
        chromosome_ends[chrom_id] = running_sum
    chromosome_ids = list(chromosome_ends.keys())
    genome_length = len(genome_str)
    return chromosome_ids, chromosome_ends, genome_str, genome_length


def generate_random_molecules(genome_file_path,
                              output_file_path=None,
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
    genome = Fasta(genome_file_path)
    chromosome_ids, chromosome_ends, genome_str, genome_length = process_genome(genome)

    for i in range(number_of_molecules):
        chromosome_idx, start, end = generate_random_molecule(molecule_size_mean,
                                                              molecule_size_standard_dev,
                                                              genome_length,
                                                              chromosome_ends,
                                                              min_molecule_size)
        chromosome = chromosome_ids[chromosome_idx]
        print(">{}:{}-{}-{}".format(i, chromosome, start, end), file=output_file)
        print(genome_str[start:end], file=output_file)

def generate_molecules_from_bed(genome_file_path,
                                bed_file_path,
                                output_file_path=None,
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
    genome_obj = Fasta(genome_file_path)
    genome = dict()
    for key in genome_obj.keys():
        genome[key] = str(genome_obj[key])
    bed_positions = list()
    for line in open(bed_file_path).readlines()[1:]:
        line = line.rstrip().split('\t')
        chr = line[0]
        start = int(line[1])
        end = int(line[2])
        gene = line[3]
        for idx in range(start, end + 1):
            bed_positions.append((chr, idx))

    num_molecules = 0
    while num_molecules < number_of_molecules:
        mol_len = math.floor(random.normalvariate(mu=molecule_size_mean, sigma=molecule_size_standard_dev))
        if mol_len < min_molecule_size:
            continue
        bed_chr, bed_pos = random.choice(bed_positions)
        if bed_pos - mol_len <= 0:
            continue
        chr_len = len(genome[bed_chr])
        if bed_pos + mol_len >= chr_len:
            continue
        mol_start = bed_pos - random.randint(0, mol_len)
        mol_end = mol_start + mol_len
        mol_str = str(genome[bed_chr][mol_start: mol_end]).upper()
        if 'N' in mol_str:
            continue
        print(">{}:{}-{}-{}".format(num_molecules, bed_chr, mol_start, mol_end), file=output_file)
        print(mol_str, file=output_file)
        num_molecules += 1

def main():
    args = parse_args()
    if args.bed:
        generate_molecules_from_bed(genome_file_path=args.reference,
                                    bed_file_path=args.bed,
                                    output_file_path=args.output_molecules,
                                    molecule_size_mean=args.molecule_size_mean,
                                    molecule_size_standard_dev=args.molecule_size_standard_dev,
                                    number_of_molecules=args.number_of_molecules,
                                    random_seed=args.random_seed,
                                    min_molecule_size=args.min_molecule_size)
    else:
        generate_random_molecules(genome_file_path=args.reference,
                                  output_file_path=args.output_molecules,
                                  molecule_size_mean=args.molecule_size_mean,
                                  molecule_size_standard_dev=args.molecule_size_standard_dev,
                                  number_of_molecules=args.number_of_molecules,
                                  random_seed=args.random_seed,
                                  min_molecule_size=args.min_molecule_size)


if __name__ == "__main__":
    main()
