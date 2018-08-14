import argparse

BLOCK_SIZE = 500

MSA_PROG_PATH = '/groups/hachgrp/projects/calib/scripts/calib_spoa/simulating/seqan_msa'
PBS_MSA_CMD = '{} {} {};'

def parse_args():
    parser = argparse.ArgumentParser(
        description="Runs MSA on server cluster")
    parser.add_argument("-f",
                        "--forward-fastq-file",
                        type=str,
                        required=True,
                        help="Input forward fastq file ")
    parser.add_argument("-r",
                        "--reverse-fastq-file",
                        type=str,
                        required=True,
                        help="Input reverse fastq file ")
    parser.add_argument("-c",
                        "--cluster-file",
                        type=str,
                        required=True,
                        help="Input cluster file ")
    parser.add_argument("-o",
                        "--output-directory",
                        type=str,
                        required=True,
                        help="Output directory ")
    args = parser.parse_args()
    return args

def create_pbs(input_file_path):
    pbs_path = '{}.pbs'.format(input_file_path)
    pbs_file = open(pbs_path, 'w+')

    print("#!/bin/bash", file=pbs_file)
    print("#SBATCH --job-name=msa", file=pbs_file)
    print("#SBATCH -c 1", file=pbs_file)
    print("#SBATCH --mem 10G", file=pbs_file)
    print("#SBATCH -t 02:59:59", file=pbs_file)
    print("#SBATCH --output={}.out".format(pbs_path), file=pbs_file)
    print("#SBATCH --error={}.err".format(pbs_path), file=pbs_file)
    print("#SBATCH --export=all", file=pbs_file)
    print("#SBATCH -p debug,express,normal,big-mem,long", file=pbs_file)
    print("{} {}".format(MSA_PROG_PATH, input_file_path), file=pbs_file)

def main():
    args = parse_args()

    cluster_to_reads = list()
    for idx, line in enumerate(open(args.cluster_file)):
        line = line.rstrip().split('\t')
        cid = int(line[0])
        rid = int(line[2])
        if not cid < len(cluster_to_reads):
            new_len = cid + 1  - len(cluster_to_reads)
            cluster_to_reads.extend([list() for _ in range(new_len)])
        cluster_to_reads[cid].append(rid)
        if idx % 1000000 == 0:
            print('Read {} cluster records; Max cluster found is {}'.format(idx, len(cluster_to_reads)))

    cluster_to_sequences = list()
    forward_file = open(args.forward_fastq_file)
    forward_mates = list()
    line = forward_file.readline()
    while len(line) > 0:
        line = forward_file.readline().rstrip().upper()
        forward_mates.append(line)
        line = forward_file.readline()
        line = forward_file.readline()
        line = forward_file.readline()
        if len(forward_mates) % 1000000 == 0:
            print('Read {} forward mates'.format(len(forward_mates)))

    reverse_file = open(args.reverse_fastq_file)
    reverse_mates = list()
    line = reverse_file.readline()
    while len(line) > 0:
        line = reverse_file.readline().rstrip().upper()
        reverse_mates.append(line)
        line = reverse_file.readline()
        line = reverse_file.readline()
        line = reverse_file.readline()
        if len(reverse_mates) % 1000000 == 0:
            print('Read {} reverse mates'.format(len(reverse_mates)))

    block_count = 100000
    for cid, reads in enumerate(cluster_to_reads):
        if cid % 100000 == 0:
            print('Processing cluster {} with {} reads'.format(cid, len(reads)))
        if cid % BLOCK_SIZE == 0:
            output_forward_file_path = '{}msa_{}_1.input'.format(args.output_directory, block_count)
            output_reverse_file_path = '{}msa_{}_2.input'.format(args.output_directory, block_count)
            create_pbs(output_forward_file_path)
            create_pbs(output_reverse_file_path)
            output_forward_file = open(output_forward_file_path, 'w+')
            output_reverse_file = open(output_reverse_file_path, 'w+')
            block_count+=1
        header = '@{}\t'.format(cid)
        for rid in cluster_to_reads[cid]:
            header += '{};'.format(rid)
        print(header, file = output_forward_file)
        print(header, file = output_reverse_file)
        for rid in reads:
            print(forward_mates[rid], file = output_forward_file)
            print(reverse_mates[rid], file = output_reverse_file)

if __name__ == "__main__":
    main()
