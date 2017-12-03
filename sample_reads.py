import numpy as np


def sample_fastq(fastq_file, choices):
    fastq = open(fastq_file)
    EOF = False
    read_is_next = False
    i = 0
    choice_index = 0
    while not EOF:
        line = fastq.readline()
        if read_is_next:
            read_is_next = False
            if choices[choice_index] == i:
                print(line)
        if line.startswith("@"):
            i += 1
            read_is_next = True
        if line == "":
            EOF = True
    fastq.close()


def main():
    sample_size = 10000000
    num_of_reads = 59998500
    fastq_file = 'amp_seq_10k_1m1.fq'
    choices = np.random.choice([i for i in range(1, num_of_reads+1)], size=sample_size, replace=False)
    sample_fastq(fastq_file, choices)


if __name__ == "__main__":
    main()