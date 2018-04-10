import sys
import argparse
import subprocess

import numpy
from pymsa.score import SumOfPairs
from pymsa.substitutionmatrix import SubstitutionMatrix

TMP_PATH = '/home/borabi/tmp/msa/'
PBS_HEADER = '#$ -N msa_{}\n'+\
                '#$ -S /bin/bash\n'+\
                '#$ -l h_rt=03:00:00\n'+\
                '#$ -l h_stack=128M\n'+\
                '#$ -l h_vmem=1G\n'+\
                '#$ -pe ncpus 1\n'+\
                '#$ -e {}.err \n'+\
                '#$ -o {}.out \n'+\
                '#$ -V\n'
BLOCK_SIZE = 500

PBS_MSA_CMD = 'make_msa {} {};'

SPS_PATH = '/home/borabi/projects/calib/calib/simulating/msa_to_sps.py'
PBS_SPS_CMD = 'python3 '+ SPS_PATH +' -f {} -r {} -o {};'

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
    parser.add_argument("-i",
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

def main():
    args = parse_args()

    forward_file = open(args.forward_fastq_file)
    forward_mates = list()
    line = forward_file.readline()
    while len(line) > 0:
        line = forward_file.readline().rstrip().upper()
        forward_mates.append(line)
        line = forward_file.readline()
        line = forward_file.readline()
        line = forward_file.readline()
    reverse_file = open(args.reverse_fastq_file)
    reverse_mates = list()
    line = reverse_file.readline()
    while len(line) > 0:
        line = reverse_file.readline().rstrip().upper()
        reverse_mates.append(line)
        line = reverse_file.readline()
        line = reverse_file.readline()
        line = reverse_file.readline()
    cluster_file = open(args.cluster_file)
    forward_clusters = list()
    reverse_clusters = list()
    print(len(forward_mates),len(reverse_mates))
    for line in cluster_file:
        if line[0] == '#':
            forward_clusters.append(list())
            reverse_clusters.append(list())
            # print(forward_clusters)
            continue
        else:
            line = line.rstrip().split('\t')
            # print(line[0], forward_mates[int(line[0])], reverse_mates[int(line[0])])
            forward_clusters[-1].append(forward_mates[int(line[1])])
            reverse_clusters[-1].append(reverse_mates[int(line[1])])
    block_count = 100001
    cluster_count = 0
    # TMP_PATH = '/home/borabi/tmp/msa/'
    # TMP_PATH += args.output_directory + '/'
    # print(TMP_PATH)
    # print(subprocess.list2cmdline(['rm','-rf',TMP_PATH]))
    # subprocess.run(['rm','-rf',TMP_PATH],shell=False)
    # subprocess.run(['mkdir',TMP_PATH],shell=False)
    block_path_prefix = TMP_PATH + args.output_directory +'_block_' +str(block_count)
    output_forward_file_path = block_path_prefix + '_1'
    output_reverse_file_path = block_path_prefix + '_2'
    output_forward_file = open(output_forward_file_path, 'w+')
    output_reverse_file = open(output_reverse_file_path, 'w+')
    for c in range(len(forward_clusters)):
        print('>{}'.format(c), file=output_forward_file)
        print('>{}'.format(c), file=output_reverse_file)
        print('\n'.join(forward_clusters[c]), file=output_forward_file)
        print('\n'.join(reverse_clusters[c]), file=output_reverse_file)
        print('<', file=output_forward_file)
        print('<', file=output_reverse_file)
        cluster_count+=1
        if (cluster_count % BLOCK_SIZE == 0) or (cluster_count == len(forward_clusters)):
            print(TMP_PATH)
            print(block_count)
            print(cluster_count)
            print(block_path_prefix)
            print(output_forward_file_path)
            print(output_reverse_file_path)

            print('Processing block {}'.format(block_count))
            pbs_path = block_path_prefix+'.sge.pbs'
            pbs_content = [PBS_HEADER.format(block_count,block_path_prefix,block_path_prefix),
                PBS_MSA_CMD.format(output_forward_file_path, output_forward_file_path+'.msa'),
                PBS_MSA_CMD.format(output_reverse_file_path, output_reverse_file_path+'.msa')]
            pbs_content = '\n'.join(pbs_content)
            pbs_file = open(pbs_path, 'w+')
            print(pbs_content, file=pbs_file)
            print(subprocess.list2cmdline(['qsub', pbs_path]))
            # subprocess.run(['qsub', pbs_path], shell=False)
            if (cluster_count < len(forward_clusters)):
                block_count +=1
                block_path_prefix = TMP_PATH + args.output_directory +'_block_' +str(block_count)
                output_forward_file_path = block_path_prefix + '_1'
                output_reverse_file_path = block_path_prefix + '_2'
                output_forward_file = open(output_forward_file_path, 'w+')
                output_reverse_file = open(output_reverse_file_path, 'w+')





if __name__ == "__main__":
    main()
