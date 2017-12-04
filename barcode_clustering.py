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
    parser.add_argument("-l", "--barcode-length", type=int, default=10,
                        help="Barcode length (default: 10)")
    parser.add_argument("-e", "--error-tolerance", type=int, default=2,
                        help="Error tolerance for barcode clustering (default: 2)")
    parser.add_argument("-m", "--minimizers-count", type=int, default=3,
                            help="Number of minimizers per read mate (default: 3)")
    parser.add_argument("-k", "--k-mer-size", type=int, default=8,
                            help="K-mer size for the minimizers (default: 8)")
    parser.add_argument("-o", "--log-file", help="Log file path.", required=True)

    args = parser.parse_args()
    return args


def template_generator(barcode_size, error_tolerance):
    template_id = 0
    for comb in combinations(range(barcode_size), barcode_size - error_tolerance):
        def template(barcode):
            return ''.join([barcode[i] for i in comb])
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
    _minimizer_count = args.minimizers_count
    _k_mer_size = args.k_mer_size

    log_file = open(args.log_file, 'w+', buffering=1)
    print('Hmmmm. Good morning?!', file=log_file)
    print('Step: Extracting barcodes...', file=log_file)
    start_time = time.time()

    # Extracting barcodes from fastq files
    barcode_lines_1 = get_barcodes(args.forward_reads, _barcode_length)
    barcode_lines_2 = get_barcodes(args.reverse_reads, _barcode_length)
    if not len(barcode_lines_1) == len(barcode_lines_2):
        print('\tYou messed up; the read files are not of the same length', file=log_file)
        sys.exit(42)
    # Getting a unique representative from each set of set of identical barcodes
    barcode_pairs_to_lines = get_barcode_pair_to_line_mapping(barcode_lines_1, barcode_lines_2)

    print('\tTotal number of unique barcode pairs:', len(barcode_pairs_to_lines), file=log_file)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)

    print('Step: Storing reverse complement of barcodes...', file=log_file)
    start_time = time.time()

    barcodes_1 = ['']*len(barcode_pairs_to_lines)
    barcodes_2 = ['']*len(barcode_pairs_to_lines)
    barcodes_rev_compl_1 = ['']*len(barcode_pairs_to_lines)
    barcodes_rev_compl_2 = ['']*len(barcode_pairs_to_lines)

    for index, barcode_pair in enumerate(barcode_pairs_to_lines.keys()):
        barcodes_1[index] = barcode_pair[0]
        barcodes_2[index] = barcode_pair[1]
        barcodes_rev_compl_1[index] = reverse_complement(barcode_pair[0])
        barcodes_rev_compl_2[index] = reverse_complement(barcode_pair[1])

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)
    print('Step: LSH of barcodes...', file=log_file)
    start_time = time.time()

    # lsh_list = [{} for _ in range(int(nCr(_barcode_length, _error_tolerance)))]
    fake_barcode = ''.join([chr(x+65) for x in range(_barcode_length)])

    barcode_graph_adjacency_sets = [{x} for x in range(len(barcode_pairs_to_lines))]

    for template, template_id in template_generator(_barcode_length, _error_tolerance):
        # TODO: Instead of lsh_list, generate adjacency set from each lsh then remove the lsh.
        lsh = dict() # lsh_list[template_id]
        print("\tTemplate {} with ID {}".format(template(fake_barcode), template_id), file=log_file)
        for barcode_num in range(len(barcodes_1)):
            barcode_1 = template(barcodes_1[barcode_num])
            barcode_1_rev = template(barcodes_rev_compl_1[barcode_num])
            barcode_2 = template(barcodes_2[barcode_num])
            barcode_2_rev = template(barcodes_rev_compl_2[barcode_num])

            new_key = barcode_1 + barcode_2
            lsh[new_key] = lsh.get(new_key, {barcode_num}).union({barcode_num})
            new_key = barcode_2_rev + barcode_1_rev
            lsh[new_key] = lsh.get(new_key, {barcode_num}).union({barcode_num})
        for val in lsh.values():
            for node in val:
                adjacent_nodes = val#.difference({node})
                barcode_graph_adjacency_sets[node].update(adjacent_nodes)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)


    print("Step: Building barcode graph from adjacency_sets...", file=log_file)
    start_time = time.time()

    barcode_graph = Graph([(node, neighbor) for node, neighbors in enumerate(barcode_graph_adjacency_sets) for neighbor in neighbors])
    print("\tGraph building from adjacency lists is completed", file=log_file)

    barcode_graph.vs['id'] = [x for x in range(len(barcode_pairs_to_lines))]
    print("\tLabeling vertices on graph is completed", file=log_file)

    barcode_graph.simplify()
    print("\tSimplifying the graph is completed", file=log_file)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)
    print("Step: Getting clusters (connected components) of the barcode graph...", file=log_file)
    start_time = time.time()

    clusters = barcode_graph.clusters()
    print('\tThere total of {} connected components'.format(len(clusters)), file=log_file)

    count = 0
    for index, connected_component in enumerate(clusters):
        vertices_count = len(connected_component)
        edges_count = 0
        for vertex in connected_component:
            edges_count += len(barcode_graph_adjacency_sets[vertex])
        edges_count -= vertices_count
        edges_count /= 2
        edges_count = int(edges_count)
        if vertices_count == 1:
            print('#{}\t{}\t{}\t{}:'.format(index, vertices_count, edges_count, 'nan'), file=log_file)
        else:
            print('#{}\t{}\t{}\t{}:'.format(index,  vertices_count, edges_count, 2*edges_count/(vertices_count * (vertices_count-1))), file=log_file)
        for barcode in connected_component:
            print(barcodes_1[barcode], barcodes_2[barcode], file=log_file)
        count += 1

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)


if __name__ == '__main__':
    main()
