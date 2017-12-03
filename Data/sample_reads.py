import random


def sample_fastq(fastq_file, choices, output):
    fastq = open(fastq_file)
    out = open(output, 'w')
    EOF = False
    read_is_next = False
    i = 0
    choice_index = 0
    while not EOF:
        line = fastq.readline()
        if read_is_next:
            read_is_next = False
            if choices[choice_index] == i:
                choice_index += 1
                print(line, file=out)
        if line.startswith("@"):
            i += 1
            read_is_next = True
        if line == "":
            EOF = True
    fastq.close()


def main():
    sample_size = 10000000
    num_of_reads = 59998500
    fastq_file1 = 'amp_seq_10k_1m1.fq'
    fastq_file2 = 'amp_seq_10k_1m2.fq'
    choices = random.sample(range(1, num_of_reads+1), k=sample_size)
    print('choices selected')
    choices.sort()
    sample_fastq(fastq_file1, choices, "amp_seq_10k_1m_sampled1.fq")
    sample_fastq(fastq_file2, choices, "amp_seq_10k_1m_sampled2.fq")


if __name__ == "__main__":
    main()
