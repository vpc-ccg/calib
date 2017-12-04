import random
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clusters barcodes by locality sensitive hashing.")
    parser.add_argument("-s", "--sample-size", type=int, help="Number of reads to sample", required=True)
    parser.add_argument("-n", "--num-of-reads", type=int, help="Total number of reads", required=True)
    parser.add_argument("-f", "--forward-reads", type=str, help="Forward read file path.", required=True)
    parser.add_argument("-r", "--reverse-reads", type=str, help="Reverse read file path.", required=True)
    return parser.parse_args()


def sample_fastq(fastq_file, choices, output):
    fastq = open(fastq_file)
    out = open(output, 'w')
    EOF = False
    i = 0
    choice_index = 0
    print_next = 0
    while not EOF:
        i += 1
        line = fastq.readline()

        if print_next > 0:
            print_next -= 1
            print(line, file=out, end='')
        if choice_index >= len(choices) and print_next == 0:
            break
        if choice_index < len(choices) and choices[choice_index]*4+1 == i:
            print_next = 3
            choice_index += 1
            print(line, file=out, end='')
        if line == "":
            EOF = True

    fastq.close()


def main():
    args = parse_args()
    sample_size = args.sample_size
    num_of_reads = args.num_of_reads
    fastq_file1 = args.forward_reads
    fastq_file2 = args.reverse_reads
    choices = random.sample(range(1, num_of_reads+1), k=sample_size)
    print('choices selected')
    print(choices)
    choices.sort()
    sample_fastq(fastq_file1, choices, args.forward_reads.split('.fq')[0]+"_sampled.fq")
    sample_fastq(fastq_file2, choices, args.reverse_reads.split('.fq')[0]+"_sampled")


if __name__ == "__main__":
    main()
