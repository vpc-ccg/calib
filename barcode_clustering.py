import sys
import time
import statistics
import numpy as np
import math
from itertools import combinations
from igraph import Graph, summary, InternalError
import argparse

_complements = {'A': 'T',
                'T': 'A',
                'G': 'C',
                'C': 'G'}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clusters barcodes by locality sensitive hashing.")
    parser.add_argument("-f", "--forward-reads", type=str, help="Forward read file path.", required=True)
    parser.add_argument("-r", "--reverse-reads", type=str, help="Reverse read file path.", required=True)
    parser.add_argument("-t", "--barcode-mini-tsv-file", type=str, help="Reverse read file path.", required=True)
    parser.add_argument("-l", "--barcode-length", type=int, default=10,
                        help="Barcode length (default: 10)")
    parser.add_argument("-e", "--error-tolerance", type=int, default=2,
                        help="Error tolerance for barcode clustering (default: 2)")
    parser.add_argument("-m", "--minimizers-count", type=int, default=3,
                            help="Number of minimizers per read mate (default: 3)")
    parser.add_argument("-x", "--minimizers-threshold", type=int, default=1,
                            help="Threshold for number of minimizers matching per read mate (default: 1)")
    # parser.add_argument("-k", "--k-mer-size", type=int, default=8,
    #                         help="K-mer size for the minimizers (default: 8)")
    parser.add_argument("-o", "--log-file", help="Log file path.", required=True)

    args = parser.parse_args()
    return args


def template_generator(barcode_size, error_tolerance):
    template_id = 0
    for comb in combinations(range(barcode_size), barcode_size - error_tolerance):
        comb = list(comb)
        def template(barcode):
            return tuple([barcode[i] for i in comb])
        yield template, template_id
        template_id += 1


def reverse_complement(sequence):
    new_seq = ""
    for char in reversed(sequence):
        new_seq += _complements.get(char, char)
    return new_seq

def nCr(n, r):
    f = math.factorial
    return int(f(n) / f(r) / f(n-r))


def ascii_to_phred(ascii):
    return 10**(-(ord(ascii) - 33)/10)


def phred_to_ascii(phred):
    if phred < 0.00000001:
        Q = 42
    else:
        Q = -10 * math.log10(phred)
    if Q > 42:
        Q = 42
    return chr(int(Q) + 33)


def consensus(sequences, qualities):
    consensus_quality = ''
    consensus = ''
    for i in range(max((len(s) for s in sequences))):
        votes = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for idx, seq in enumerate(sequences):
            if len(seq) <= i:
                continue
            votes[seq[i]] += math.log10(1-ascii_to_phred(qualities[idx][i]))
            for b in votes:
                if b == seq[i]:
                    continue
                votes[b] += math.log10(ascii_to_phred(qualities[idx][i]))
        for v in votes:
            votes[v] = math.pow(10, votes[v])
        norm_fact = sum(votes.values())
        for v in votes:
            votes[v] = votes[v]/norm_fact
        #print(votes)
        consensus += max(votes.keys(), key=(lambda k: votes[k]))
        consensus_quality += phred_to_ascii(1-votes[consensus[-1]])
    return consensus, consensus_quality

def fastq_lines(fastq_path):
    fastq_file = open(fastq_path)
    fastq_lines = fastq_file.readlines()
    read_count = int(len(fastq_lines)/4)
    fastq_lines_tuples = [0]*read_count
    for idx in range(read_count):
        fastq_line_idx = idx*4
        fastq_lines_tuples[idx] = (fastq_lines[fastq_line_idx].rstrip(), fastq_lines[fastq_line_idx+1].rstrip(), fastq_lines[fastq_line_idx+3].rstrip())
    return fastq_lines_tuples

dna_encoding = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
def dna_to_index(DNA):
    index = 0
    power = 1
    for letter in DNA:
        index = index + dna_encoding[letter]*power
        power = power*4
    return index


def main():
    # Parsing args
    args = parse_args()
    print(args)
    _barcode_length = args.barcode_length
    _error_tolerance = args.error_tolerance
    _barcode_mini_tsv = args.barcode_mini_tsv_file
    _minimizer_count = args.minimizers_count
    _minimizers_threshold = args.minimizers_threshold
    # _k_mer_size = args.k_mer_size

    log_file = open(args.log_file, 'w+', buffering=1)
    print('Hmmmm. Good morning?!', file=log_file)
    print('Step: Extracting barcodes...', file=log_file)
    start_time = time.time()

    node_to_reads_dict = dict()

    _barcode_mini_tsv_file = open(_barcode_mini_tsv)
    for idx, line in enumerate(_barcode_mini_tsv_file):
        node_key = tuple(line.rstrip().split('\t'))
        if node_key in node_to_reads_dict:
            node_to_reads_dict[node_key].append(idx)
        else:
            node_to_reads_dict[node_key] = [idx]

    node_count = len(node_to_reads_dict)
    node_to_reads = np.zeros(shape=node_count, dtype=object)
    node_to_barcode = np.empty(shape=node_count, dtype='S{}'.format(_barcode_length))
    node_to_mini_1 = np.zeros(shape=(node_count, _minimizer_count), dtype=np.int64)
    node_to_mini_2 = np.zeros(shape=(node_count, _minimizer_count), dtype=np.int64)

    for idx, key in enumerate(node_to_reads_dict):
        node_to_reads[idx] = node_to_reads_dict[key]
        node_to_barcode[idx] = key[0]
        node_to_mini_1[idx] = key[1:1+_minimizer_count]
        node_to_mini_2[idx] = key[1+_minimizer_count:1+_minimizer_count+_minimizer_count]
    del(idx)
    del(key)


    # for idx in range(node_count):
    #     print(node_to_reads[idx])
    #     print(node_to_barcode_1[idx], node_to_barcode_2[idx], node_to_mini_1[idx], node_to_mini_2[idx], sep='\t')

    print('\tTotal number of nodes:', node_count, file=log_file)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)

    print('Step: LSH of barcodes...', file=log_file)
    start_time = time.time()

    adjacency_sets = [{x} for x in range(node_count)]

    for template, template_id in template_generator(_barcode_length, _error_tolerance):
        lsh = dict()
        for idx in range(node_count):
            templated_barcode = template(node_to_barcode[idx])
            if 'N' in templated_barcode:
                continue
            if templated_barcode in lsh:
                lsh[templated_barcode].append(idx)
            else:
                lsh[templated_barcode] = [idx]
        for key in lsh:
            lsh[key] = set(lsh[key])
            adjacent_nodes = lsh[key]
            for node in adjacent_nodes:
                adjacency_sets[node].update(adjacent_nodes)
        del(lsh)
    # for idx, neighbors in enumerate(adjacency_sets):
    #     print(idx, neighbors, sep='\t')
    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)

    print("Step: Removing edges due to mini\'s...", file=log_file)
    start_time = time.time()

    for node, neighbors in enumerate(adjacency_sets):
        neighbors.remove(node)
        for neighbor in list(neighbors):
            # print(node_to_mini_1[node], node_to_mini_1[neighbor], sep='\t')
            # print(node_to_mini_2[node], node_to_mini_2[neighbor], sep='\t')
            # print([node_to_mini_1[node][i]==node_to_mini_1[neighbor][i] for i in range(_minimizer_count)])
            # print([node_to_mini_2[node][i]==node_to_mini_2[neighbor][i] for i in range(_minimizer_count)])
            matching_minimizers_1 = sum([node_to_mini_1[node][i]==node_to_mini_1[neighbor][i] for i in range(_minimizer_count)])
            matching_minimizers_2 = sum([node_to_mini_2[node][i]==node_to_mini_2[neighbor][i] for i in range(_minimizer_count)])
            if (matching_minimizers_1 < _minimizers_threshold) or (matching_minimizers_2 < _minimizers_threshold):
                neighbors.remove(neighbor)
                adjacency_sets[neighbor].remove(node)
    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)


    print("Step: Building graph from adjacency_sets...", file=log_file)
    start_time = time.time()

    graph = Graph([(node, neighbor) for node, neighbors in enumerate(adjacency_sets) for neighbor in neighbors])
    print("\tGraph building from adjacency lists is completed", file=log_file)

    graph.vs['id'] = [x for x in range(node_count)  ]
    print("\tLabeling vertices on graph is completed", file=log_file)

    graph.simplify()
    print("\tSimplifying the graph is completed", file=log_file)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)
    print("Step: Getting clusters (connected components) of the barcode graph...", file=log_file)
    start_time = time.time()

    clusters = graph.clusters()
    print('\tThere total of {} connected components'.format(len(clusters)), file=log_file)

    fastq_1 = fastq_lines(args.forward_reads)
    fastq_2 = fastq_lines(args.reverse_reads)
    out_fastq_1 = open(args.forward_reads + '.corrected.fastq', 'w+')
    out_fastq_2 = open(args.reverse_reads + '.corrected.fastq', 'w+')

    if len(fastq_1) != len(fastq_2):
        print('You\'ve screwed up; fastq files don\'t match in line count')
        exit()

    for index, connected_component in enumerate(clusters):
        nodes_count = len(connected_component)
        edges_count = 0
        read_count = 0
        # mate_1_tuples = []
        # mate_2_tuples = []
        for node in connected_component:
            edges_count += len(adjacency_sets[node])
            read_count += len(node_to_reads[node])
            # for read in node_to_reads[node]:
            #     mate_1_tuples.append(fastq_1[read])
            #     mate_2_tuples.append(fastq_2[read])
        # edges_count -= nodes_count
        edges_count /= 2
        edges_count = int(edges_count)
        if nodes_count == 1:
            print('#{}\t{}\t{}\t{}\t{}:'.format(index, nodes_count, edges_count, read_count,'nan'), file=log_file)
        else:
            print('#{}\t{}\t{}\t{}\t{}:'.format(index,  nodes_count, edges_count, read_count, 2*edges_count/(nodes_count * (nodes_count-1))), file=log_file)
        for node in connected_component:
            print(node_to_barcode[node], node_to_mini_1[node], node_to_mini_2[node], file=log_file)
        # seq_1, qual_1 = consensus(sequences=[i[1] for i in mate_1_tuples], qualities=[i[2] for i in mate_1_tuples])
        # seq_2, qual_2 = consensus(sequences=[i[1] for i in mate_2_tuples], qualities=[i[2] for i in mate_2_tuples])
        # for read in range(read_count):
        #     print(fastq_1[read][0], file=out_fastq_1)
        #     print(fastq_1[read][1], file=out_fastq_1)
        #     print('+', file=out_fastq_1)
        #     print(fastq_1[read][2], file=out_fastq_1)
        #
        #     print(fastq_2[read][0], file=out_fastq_2)
        #     print(fastq_2[read][1], file=out_fastq_2)
        #     print('+', file=out_fastq_2)
        #     print(fastq_2[read][2], file=out_fastq_2)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)


if __name__ == '__main__':
    main()
