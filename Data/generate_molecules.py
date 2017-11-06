from pyfaidx import Fasta

import random
import math

molecule_size_mu = 150
molecule_size_sig = 20
number_of_molecules = 100
random_seed = 42
output_file_path = "molecules.fa"
genome = Fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")


random.seed(random_seed)

genome_length = len(genome[0])
with open(output_file_path, 'w+') as output_file:
    for i in range(number_of_molecules):
        start = random.randrange(0, genome_length)
        end = min([start + math.floor(random.normalvariate(mu=molecule_size_mu, sigma=molecule_size_sig)), genome_length])
        print('>'+str(i)+':'+str(start)+'-'+str(end), file=output_file)
        print(genome[0][start:end], file=output_file)
