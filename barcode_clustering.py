import sys
import time
import numpy as np
import math
from itertools import combinations
from igraph import Graph, summary, InternalError
import argparse


_complements = {'A': 'T',
                'T': 'A',
                'G': 'C',
                'C': 'G'}


def reverse_complement(sequence):
    new_seq = ""
    for char in reversed(sequence):
        new_seq += _complements.get(char, char)
    return new_seq


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clusters barcodes by locality sensitive hashing.")
    parser.add_argument("-f",
                        "--forward-reads",
                        type=str,
                        required=True,
                        help="Forward read file path.")
    parser.add_argument("-r",
                        "--reverse-reads",
                        type=str,
                        required=True,
                        help="Reverse read file path.")
    parser.add_argument("-t",
                        "--barcode-mini-tsv-file",
                        type=str,
                        required=True,
                        help="Barcode and minimizers tab seperated values file.")
    parser.add_argument("-e",
                        "--barcode-error-tolerance",
                        type=int,
                        default=2,
                        help="Error tolerance for barcode clustering (default: 2)")
    parser.add_argument("-m",
                        "--minimizers-threshold",
                        type=int,
                        default=1,
                        help="Threshold for number of minimizers matching per read mate (default: 1)")
    parser.add_argument("-q",
                        "--ratio",
                         type=float,
                         default=1.0,
                         help="The ratio of nodes to use per template (default: 1.0)")
    parser.add_argument("-s",
                        "--random-seed",
                        type=int,
                        default=42,
                        help="NumPy random seed (default: 42)")
    parser.add_argument("-o",
                        "--output-prefix",
                        required=True,
                        help="Output prefix for log, clusters, and collapesd (corrected) fastq files.")
    parser.add_argument("-c",
                        "--correct",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help="Output corrected fastq files")
    args = parser.parse_args()
    return args

# Given a barcode size and error tolerance, generates C(barcode size, error tolerance) functions.
# Each function takes in a barcode, and returns the same barcode with error-tolerance characters dropped
# This allows for comparing barcodes on all possible templates to match barcodes that are similar within
# the specified error-tolerance.
def template_generator(barcode_size, error_tolerance):
    l = list(enumerate(combinations(range(barcode_size), barcode_size - error_tolerance)))
    np.random.shuffle(l)
    for template_id, comb in l:
        comb = list(comb)
        def template(barcode):
            return tuple([barcode[i] for i in comb])
        yield template, template_id


def nCr(n, r):
    f = math.factorial
    return int(f(n) / f(r) / f(n-r))

# Converts FASTQ ASCII quality to error probablity
def ascii_to_phred(ascii):
    return 10**(-(ord(ascii) - 33)/10)

# Converts error probablity to FASTQ ASCII quality
def phred_to_ascii(phred):
    if phred < 0.00000001:
        Q = 42
    else:
        Q = -10 * math.log10(phred)
    if Q > 42:
        Q = 42
    return chr(int(Q) + 33)

# Given a set of sequences and their quality scores, produces the sequence with the most
# maximum likelihood with its quality score its likelihood
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
        consensus += max(votes.keys(), key=(lambda k: votes[k]))
        consensus_quality += phred_to_ascii(1-votes[consensus[-1]])
    return consensus, consensus_quality


# Loads a FASTQ file into a list of 3-tuples (read name, sequence, quality)
def fastq_lines(fastq_path):
    fastq_file = open(fastq_path)
    fastq_lines = fastq_file.readlines()
    read_count = int(len(fastq_lines)/4)
    fastq_lines_tuples = [0]*read_count
    for idx in range(read_count):
        fastq_line_idx = idx*4
        fastq_lines_tuples[idx] = (fastq_lines[fastq_line_idx].rstrip(), fastq_lines[fastq_line_idx+1].rstrip(), fastq_lines[fastq_line_idx+3].rstrip())
    return fastq_lines_tuples


def main():
    # Parsing args
    args = parse_args()
    _barcode_mini_tsv = args.barcode_mini_tsv_file
    _error_tolerance = args.barcode_error_tolerance
    _minimizers_threshold = args.minimizers_threshold
    _random_seed = args.random_seed
    _output_prefix = args.output_prefix

    # The log file has buffering of one line, to ensure that it is updated frequently
    log_file = open("{}.log".format(_output_prefix), 'w+', buffering=1)

    # Random seed for reproducibility
    np.random.seed(_random_seed)

    print('Hmmmm. Good morning?!', file=log_file)
    print('Hmmmm. Good morning?!')
    print('Well, it is morning somewhere...', file=log_file)
    print('Well, it is morning somewhere...')
    print('So good morning it is!', file=log_file)
    print('So good morning it is!')
    print('Step: Extracting barcodes...', file=log_file)
    print('Step: Extracting barcodes...')
    start_time = time.time()

    # Maps each node to set of reads that have identical barcodes and minimizers
    node_to_reads_dict = dict()
    barcode_mini_tsv_file = open(_barcode_mini_tsv)
    for idx, line in enumerate(barcode_mini_tsv_file):
        node_key = tuple(line.rstrip().split('\t'))
        if node_key in node_to_reads_dict:
            node_to_reads_dict[node_key].append(idx)
        else:
            node_to_reads_dict[node_key] = [idx]
    _barcode_length = len(node_key[0])
    _minimizer_count = (len(node_key) - 1 ) // 2

    # Checking that each mate as the same number of minimizers
    if (len(node_key) -1) % 2 != 0:
        print('Something is wrong with your TSV file', file=log_file)
        print('Something is wrong with your TSV file')
        exit()

    node_count = len(node_to_reads_dict)
    print('\tTotal number of nodes:', node_count, file=log_file)
    print('\tTotal number of nodes:', node_count)

    # Turning the node_to_reads_dict dictionary into 4 lists, then deleting the dictionary
    node_to_reads = np.zeros(shape=node_count, dtype=object)
    node_to_barcode = np.empty(shape=node_count, dtype='S{}'.format(_barcode_length))
    node_to_mini_1 = np.zeros(shape=(node_count, _minimizer_count), dtype=np.int64)
    node_to_mini_2 = np.zeros(shape=(node_count, _minimizer_count), dtype=np.int64)
    for idx, key in enumerate(node_to_reads_dict):
        node_to_reads[idx] = node_to_reads_dict[key]
        node_to_barcode[idx] = key[0]
        node_to_mini_1[idx] = key[1:1+_minimizer_count]
        node_to_mini_2[idx] = key[1+_minimizer_count:1+_minimizer_count+_minimizer_count]
    del(node_to_reads_dict)
    del(idx)
    del(key)


    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)
    print('\tLast step took {} seconds'.format(finish_time - start_time))

    print('Step: LSH of barcodes...', file=log_file)
    print('Step: LSH of barcodes...')
    start_time = time.time()

    # We now want to build the graph. We start by adding edges implied by barcode similarity
    adjacency_sets = [{x} for x in range(node_count)]
    # We iterate over each possible template on which barcodes may match on
    for tempalte_order, (template, template_id) in enumerate(template_generator(_barcode_length, _error_tolerance)):
        lsh = dict()
        # We choose a random subset of the reads for this tempalte to process.
        # The size of the subset is defined by the user.
        for idx in np.random.choice(node_count, size= max(1, int(node_count*args.ratio)) ,replace=False):
            templated_barcode = template(node_to_barcode[idx])
            # No matches on N's are allowed
            if 'N' in templated_barcode:
                continue
            if templated_barcode in lsh:
                lsh[templated_barcode].append(idx)
            else:
                lsh[templated_barcode] = [idx]
        # Report how many new edges this template has added
        count_new_edges = 0
        # All nodes that sharea bucket should be neighbors of one another
        for key in lsh:
            lsh[key] = set(lsh[key])
            adjacent_nodes = lsh[key]
            for node in adjacent_nodes:
                count_new_edges += len(adjacent_nodes.difference(adjacency_sets[node]))
                adjacency_sets[node].update(adjacent_nodes)
        print('\t{}: template {} added {} new edges'.format(tempalte_order, template_id, count_new_edges), file=log_file)
        print('\t{}: template {} added {} new edges'.format(tempalte_order, template_id, count_new_edges))
        del(lsh)
    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)
    print('\tLast step took {} seconds'.format(finish_time - start_time))

    print("Step: Removing edges due to mini\'s...", file=log_file)
    print("Step: Removing edges due to mini\'s...")
    start_time = time.time()

    # Now remove any neighbors that don't have enough matching minimizers on the left or the right
    for node, neighbors in enumerate(adjacency_sets):
        neighbors.remove(node)
        for neighbor in list(neighbors):
            matching_minimizers_1 = sum([node_to_mini_1[node][i]==node_to_mini_1[neighbor][i] for i in range(_minimizer_count)])
            matching_minimizers_2 = sum([node_to_mini_2[node][i]==node_to_mini_2[neighbor][i] for i in range(_minimizer_count)])
            if (matching_minimizers_1 < _minimizers_threshold) or (matching_minimizers_2 < _minimizers_threshold):
                neighbors.remove(neighbor)
                adjacency_sets[neighbor].remove(node)
    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)
    print('\tLast step took {} seconds'.format(finish_time - start_time))


    print("Step: Building graph from adjacency_sets...", file=log_file)
    print("Step: Building graph from adjacency_sets...")
    start_time = time.time()

    # From the adjacency sets, create an igraph graph
    graph = Graph([(node, neighbor) for node, neighbors in enumerate(adjacency_sets) for neighbor in neighbors])
    print("\tGraph building from adjacency lists is completed", file=log_file)
    print("\tGraph building from adjacency lists is completed")

    # Labeling the nodes
    graph.vs['id'] = [x for x in range(node_count)  ]
    print("\tLabeling vertices on graph is completed", file=log_file)
    print("\tLabeling vertices on graph is completed")

    # Remove any loops
    graph.simplify()
    print("\tSimplifying the graph is completed", file=log_file)
    print("\tSimplifying the graph is completed")

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)
    print('\tLast step took {} seconds'.format(finish_time - start_time))
    print("Step: Getting clusters (connected components) of the barcode graph...", file=log_file)
    print("Step: Getting clusters (connected components) of the barcode graph...")
    start_time = time.time()

    # Get connected components of the graph
    clusters = graph.clusters()
    print('\tThere total of {} connected components'.format(len(clusters)), file=log_file)
    print('\tThere total of {} connected components'.format(len(clusters)))

    fastq_1 = fastq_lines(args.forward_reads)
    fastq_2 = fastq_lines(args.reverse_reads)

    if len(fastq_1) != len(fastq_2):
        print('You\'ve screwed up -_-; fastq files line counts don\'t match', file=log_file)
        print('You\'ve screwed up -_-; fastq files line counts don\'t match')
        exit()

    # Process each cluster by outputting it with its statistics in .clusters file
    clusters_file = open(_output_prefix + '.clusters', 'w+')
    for index, connected_component in enumerate(clusters):
        nodes_count = len(connected_component)
        edges_count = 0
        read_count = 0
        reads = []
        for node in connected_component:
            edges_count += len(adjacency_sets[node])
            reads.extend(node_to_reads[node])
        # Don't forget to count A-B and B-A edges as one edge (this is an undirected graph)
        edges_count = edges_count // 2
        read_count = len(reads)
        # Report #cluster_id node_count edge_count read_count density for each cluster
        if nodes_count == 1:
            print('#{}\t{}\t{}\t{}\t{}:'.format(index, nodes_count, edges_count, read_count,'nan'), file=clusters_file)
        else:
            print('#{}\t{}\t{}\t{}\t{}:'.format(index,  nodes_count, edges_count, read_count, 2*edges_count/(nodes_count * (nodes_count-1))), file=clusters_file)
        # Print the reads in each cluster
        for read in reads:
            print(fastq_1[read][0], fastq_1[read][1], fastq_1[read][2], fastq_2[read][0], fastq_2[read][1], fastq_2[read][2], sep='\t', file=clusters_file)
    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)

    if not args.correct:
        exit()

    print("Step: Correcting and collapsing clusters...", file=log_file)
    start_time = time.time()

    out_fastq_1 = open(_output_prefix + '.collapsed_1.fq', 'w+')
    out_fastq_2 = open(_output_prefix + '.collapsed_2.fq', 'w+')
    # For each cluster, output one concesus sequence with its quality score
    for index, connected_component in enumerate(clusters):
        sequences_1 = []
        qualities_1 = []
        sequences_2 = []
        qualities_2 = []
        for node in connected_component:
            for read in node_to_reads[node]:
                sequences_1.append(fastq_1[read][1])
                qualities_1.append(fastq_1[read][2])
                sequences_2.append(fastq_2[read][1])
                qualities_2.append(fastq_2[read][2])
        seq_1, qual_1 = consensus(sequences=sequences_1, qualities=qualities_1)
        print('@cluster_id_{}/1\n{}\n{}\n{}'.format(index, seq_1, '+', qual_1), file=out_fastq_1)
        seq_2, qual_2 = consensus(sequences=sequences_2, qualities=qualities_2)
        print('@cluster_id_{}/2\n{}\n{}\n{}'.format(index, seq_2, '+', qual_2), file=out_fastq_2)

    finish_time = time.time()
    print('\tLast step took {} seconds'.format(finish_time - start_time), file=log_file)


if __name__ == '__main__':
    main()
