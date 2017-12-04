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
    parser.add_argument("-k", "--k-mer-size", type=int, default=8,
                            help="K-mer size for the minimizers (default: 8)")
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


def get_barcodes(fastq_file, barcode_length):
    fastq = open(fastq_file)
    barcodes = []
    EOF = False
    read_is_next = False
    while not EOF:
        line = fastq.readline()
        if read_is_next:
            barcodes.append(line[0:barcode_length])
            read_is_next = False
        if line.startswith("@"):
            read_is_next = True
        if line == "":
            EOF = True
    fastq.close()
    return barcodes


def get_barcode_pair_to_line_mapping(barcode_lines_1, barcode_lines_2):
    barcode_pair_to_line_dict = dict()
    for line in range(len(barcode_lines_1)):
        barcode_pair = (barcode_lines_1[line], barcode_lines_2[line])
        if barcode_pair in barcode_pair_to_line_dict:
            barcode_pair_to_line_dict[barcode_pair].add(line)
        else:
            barcode_pair_to_line_dict[barcode_pair] = {line}
    return barcode_pair_to_line_dict


def get_lsh(barcodes, error_tolerance):
    barcode_length = len(barcodes[0])
    # From AAAAAAAA, AAAAAAAT, ..., TTTTTTTT, ..., CCCCCCCC, ..., GGGGGGGG
    # lsh = np.array([set()]*(4**(barcode_length-error_tolerance)))
    # TODO: Implement this without hashing using the above method
    lsh = {}
    for template, template_id in template_generator(barcode_length, error_tolerance):
        for index in range(len(barcodes)):
            barcode_word = template_id + template(barcodes[index])
            lsh[barcode_word] = lsh.get(barcode_word, {index}).union({index})
    return lsh


def get_adjacency_set(lsh, num_barcodes):
    # TODO: Once get_lsh is implemented without hashing, change function to use list instead of dictionary
    adjacency_sets = np.array([set()] * num_barcodes)
    for adjacent_elements in lsh.values():
        for node in adjacent_elements:
            adjacency_sets[node].update(lsh.values)
    return adjacency_sets


def nCr(n, r):
    f = math.factorial
    return int(f(n) / f(r) / f(n-r))


def ascii_to_phred(ascii):
    return 10**(-(ord(ascii) - 33)/10)


def phred_to_ascii(phred):
    Q = -10 * math.log10(phred)
    return chr(int(Q) + 33)


def consensus(sequences, qualities):
    # TODO: Assuming sequences of same length
    # TODO: Returns numpy chararray. Might have to process for printing.
    sequence_length = len(sequences[0])
    consensus_seq = np.chararray(shape=sequence_length)
    for i in range(sequence_length):
        probabilities_of_error = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        for sequence, quality in zip(sequences, qualities):
            nt = sequence[i]
            ascii_qual = quality[i]
            # P(error) = P(error in each column)
            probabilities_of_error[nt] *= ascii_to_phred(ascii_qual)
        consensus_seq[i] = min(probabilities_of_error.keys(), key=(lambda k: probabilities_of_error[k]))
    return consensus_seq


def main():
    # Parsing args
    args = parse_args()
    _barcode_length = args.barcode_length
    _error_tolerance = args.error_tolerance
    _barcode_mini_tsv = args.barcode_mini_tsv_file
    _minimizer_count = args.minimizers_count
    _minimizers_threshold = args.minimizers_threshold
    _k_mer_size = args.k_mer_size

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
    node_to_barcode_1 = np.empty(shape=node_count, dtype='S{}'.format(_barcode_length))
    node_to_barcode_2 = np.empty(shape=node_count, dtype='S{}'.format(_barcode_length))
    node_to_mini_1 = np.zeros(shape=(node_count, _minimizer_count), dtype=np.int64)
    node_to_mini_2 = np.zeros(shape=(node_count, _minimizer_count), dtype=np.int64)

    for idx, key in enumerate(node_to_reads_dict):
        node_to_reads[idx] = node_to_reads_dict[key]
        node_to_barcode_1[idx] = key[0]
        node_to_barcode_2[idx] = key[1]
        node_to_mini_1[idx] = key[2:2+_minimizer_count]
        node_to_mini_2[idx] = key[2+_minimizer_count:2+_minimizer_count+_minimizer_count]
    # for idx in range(node_count):
    #     print(node_to_reads[idx])
    #     print(node_to_barcode_1[idx], node_to_barcode_2[idx], node_to_mini_1[idx], node_to_mini_2[idx], sep='\t')

    print('\tTotal number of nodes:', node_count, file=log_file)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)

    print('Step: LSH of barcodes...', file=log_file)
    start_time = time.time()

    fake_barcode = ''.join([chr(x+65) for x in range(_barcode_length)])

    adjacency_sets = [{x} for x in range(node_count)]

    for template, template_id in template_generator(_barcode_length, _error_tolerance):
        # TODO: Instead of lsh_list, generate adjacency set from each lsh then remove the lsh.
        lsh = dict() # lsh_list[template_id]
        # print("\tTemplate {} with ID {}".format(template(fake_barcode), template_id), file=log_file)
        for idx in range(node_count):
            barcode_1 = template(node_to_barcode_1[idx])
            barcode_2 = template(node_to_barcode_2[idx])
            key = barcode_1 + barcode_2
            if key in lsh:
                lsh[key].append(idx)
            else:
                lsh[key] = [idx]
        for adjacent_nodes in lsh.values():
            for node in adjacent_nodes:
                adjacency_sets[node].update(set(adjacent_nodes))
    # for idx, neighbors in enumerate(adjacency_sets):
    #     print(idx, neighbors, sep='\t')
    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)

    print("Step: Removing edges due to mini\'s...", file=log_file)
    start_time = time.time()

    for node, neighbors in enumerate(adjacency_sets):
        for neighbor in list(neighbors):
            # print(node_to_mini_1[node], node_to_mini_1[neighbor], sep='\t')
            # print(node_to_mini_2[node], node_to_mini_2[neighbor], sep='\t')
            # print([node_to_mini_1[node][i]==node_to_mini_1[neighbor][i] for i in range(_minimizer_count)])
            # print([node_to_mini_2[node][i]==node_to_mini_2[neighbor][i] for i in range(_minimizer_count)])
            if (sum([node_to_mini_1[node][i]==node_to_mini_1[neighbor][i] for i in range(_minimizer_count)]) < _minimizers_threshold) or \
            ((sum([node_to_mini_2[node][i]==node_to_mini_2[neighbor][i] for i in range(_minimizer_count)]) < _minimizers_threshold)):
                neighbors.remove(neighbor)
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

    
    for index, connected_component in enumerate(clusters):
        nodes_count = len(connected_component)
        edges_count = 0
        read_count = 0
        for node in connected_component:
            edges_count += len(adjacency_sets[node])
            read_count += len(node_to_reads[node])
        edges_count -= nodes_count
        edges_count /= 2
        edges_count = int(edges_count)
        if nodes_count == 1:
            print('#{}\t{}\t{}\t{}\t{}:'.format(index, nodes_count, edges_count, read_count,'nan'), file=log_file)
        else:
            print('#{}\t{}\t{}\t{}\t{}:'.format(index,  nodes_count, edges_count, read_count, 2*edges_count/(nodes_count * (nodes_count-1))), file=log_file)
        for node in connected_component:
            print(node_to_barcode_1[node], node_to_barcode_2[node], node_to_mini_1[node], node_to_mini_2[node], file=log_file)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)


if __name__ == '__main__':
    main()
