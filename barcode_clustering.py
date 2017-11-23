import sys
import time
import statistics
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
        barcode_pair = (barcode_lines_1[line],barcode_lines_2[line])
        if barcode_pair in barcode_pair_to_line_dict:
            barcode_pair_to_line_dict[barcode_pair].add(line)
        else:
            barcode_pair_to_line_dict[barcode_pair] = {line}
    return barcode_pair_to_line_dict


def breadth_first_traversal(graph, node, cluster):
    cluster.add(node)
    nodes_to_traverse = graph[node] - cluster
    cluster.update(graph[node])
    for node in nodes_to_traverse:
        cluster.update(breadth_first_traversal(graph, node, cluster))
    return cluster


def generate_clusters_by_bfs(graph):
    # graph is a dict, clusters is an array of sets
    clusters = []
    for node in graph.keys():
        node_already_in_clusters = False
        for cluster in clusters:
            if node in cluster:
                node_already_in_clusters = True
                cluster.update(graph[node])
                break
        if not node_already_in_clusters:
            clusters.append(breadth_first_traversal(graph, node, set()))
    return clusters

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


def main():
    args = parse_args()

    _barcode_length = args.barcode_length
    _error_tolerance = args.error_tolerance

    log_file = open(args.log_file, 'a')
    print('Hmmmm. Good morning?!', file=log_file)
    print('Step: Extracting barcodes...', file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    start_time = time.time()

    barcode_lines_1 = get_barcodes(args.forward_reads, _barcode_length)
    barcode_lines_2 = get_barcodes(args.reverse_reads, _barcode_length)
    if not len(barcode_lines_1) == len(barcode_lines_2):
        log_file = open(args.log_file, 'a')
        print('You messed up; the read files are not of the same length', file=log_file)
        if not log_file == sys.stdout:
            log_file.close()
        sys.exit(42)
    barcode_pairs_to_lines = get_barcode_pair_to_line_mapping(barcode_lines_1, barcode_lines_2)

    log_file = open(args.log_file, 'a')
    print('Total number of unique barcode pairs:', len(barcode_pairs_to_lines), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()

    finish_time = time.time()
    log_file = open(args.log_file, 'a')
    print('Last step took {} seconds'.format(finish_time - start_time), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()



    # for p in barcode_pairs_to_lines:
    #     print(p, len(barcode_pairs_to_lines[p]), barcode_pairs_to_lines[p])


    log_file = open(args.log_file, 'a')
    print('Step: Storing reverse complement of barcodes...', file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
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
    log_file = open(args.log_file, 'a')
    print('Last step took {} seconds'.format(finish_time - start_time), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    log_file = open(args.log_file, 'a')
    print('Step: LSH of barcodes...', file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    start_time = time.time()

    lsh_list = [{} for _ in range(int(nCr(_barcode_length, _error_tolerance)))]
    fake_barcode = ''.join([chr(x+65) for x in range(_barcode_length)])


    for template, template_id in template_generator(_barcode_length, _error_tolerance):
        log_file = open(args.log_file, 'a')
        lsh = lsh_list[template_id]
        print("Template {} with ID {}".format(template(fake_barcode), template_id), file=log_file)
        if not log_file == sys.stdout:
            log_file.close()
        for barcode_num in range(len(barcodes_1)):
            barcode_1 = template(barcodes_1[barcode_num])
            barcode_1_rev = template(barcodes_rev_compl_1[barcode_num])
            barcode_2 = template(barcodes_2[barcode_num])
            barcode_2_rev = template(barcodes_rev_compl_2[barcode_num])

            new_key = barcode_1 + barcode_2
            lsh[new_key] = lsh.get(new_key, {barcode_num}).union({barcode_num})
            new_key = barcode_2_rev + barcode_1_rev
            lsh[new_key] = lsh.get(new_key, {barcode_num}).union({barcode_num})

    log_file = open(args.log_file, 'a')
    print('There are {} buckets with values in the LSH dictionaries'.format(sum(( len(lsh.keys()) for lsh in lsh_list))), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    finish_time = time.time()
    log_file = open(args.log_file, 'a')
    print('Last step took {} seconds'.format(finish_time - start_time), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    log_file = open(sys.argv[3] + '.supp', 'w+')
    print('Step: Counting the size of each set in LSH dictionaries...', file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    start_time = time.time()


    count = 0
    log_file = open(sys.argv[3] + '.supp', 'a')
    for lsh in lsh_list:
        for _set in lsh.values():
            print(len(_set), _set, sep='\t', file=log_file)
            count += 1
            if count == 100000:
                count = 0
                if not log_file == sys.stdout:
                    log_file.close()
                log_file = open(sys.argv[3] + '.supp', 'a')
    if not log_file == sys.stdout:
        log_file.close()
    finish_time = time.time()
    log_file = open(sys.argv[3] + '.supp', 'a')
    print('Last step took {} seconds'.format(finish_time - start_time), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    log_file = open(args.log_file, 'a')
    print("Step: Building barcode graph adjacency sets...", file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    start_time = time.time()

    barcode_graph_adjacency_sets = dict()
    count = 0
    for lsh in lsh_list:
        for val in lsh.values():
            count += 1
            if count % 100000 == 0:
                log_file = open(args.log_file, 'a')
                print('Count', count, file=log_file)
            if not log_file == sys.stdout:
                log_file.close()
            for node in val:
                adjacent_nodes = val#.difference({node})
                if node in barcode_graph_adjacency_sets:
                    barcode_graph_adjacency_sets[node].update(adjacent_nodes)
                else:
                    barcode_graph_adjacency_sets[node] = adjacent_nodes

    finish_time = time.time()
    log_file = open(args.log_file, 'a')
    print('Last step took {} seconds'.format(finish_time - start_time), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    log_file = open(args.log_file, 'a')
    print("Step: Building barcode graph from adjacency_sets...", file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    start_time = time.time()


    barcode_graph = Graph()
    barcode_graph.add_vertices(len(barcode_pairs_to_lines))
    barcode_graph.vs['id'] = [x for x in range(len(barcode_pairs_to_lines))]

    # barcode_pair_length = _barcode_length*2 - _error_tolerance*2
    count = 0
    for node in barcode_graph_adjacency_sets:
        barcode_graph.add_edges([(node, neighbor) for neighbor in barcode_graph_adjacency_sets[node]])
        count += 1
        if count % 100000 == 0:
            log_file = open(args.log_file, 'a')
            print('Count', count, file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    barcode_graph.simplify()
    finish_time = time.time()
    log_file = open(args.log_file, 'a')
    print('Last step took {} seconds'.format(finish_time - start_time), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    log_file = open(args.log_file, 'a')
    print("Step: Getting clusters (connected components) of the barcode graph...", file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    start_time = time.time()

    cc_count = 0
    for connected_component in barcode_graph.clusters().subgraphs():
        cc_count += 1
        log_file = open(args.log_file, 'a')
        print('===\nDensity: {} and vertices:'.format(connected_component.density()), file=log_file)
        for barcode in connected_component.vs['id']:
            print(barcodes_1[barcode], barcodes_2[barcode], file=log_file)
        log_file.close()
    log_file = open(args.log_file, 'a')
    print('There total of {} connected components'.format(cc_count), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    finish_time = time.time()
    log_file = open(args.log_file, 'a')
    print('Last step took {} seconds'.format(finish_time - start_time), file=log_file)
    if not log_file == sys.stdout:
        log_file.close()
    #print([len(cluster) for cluster in clusters])

    # print('len of barcode_graph {}'.format(len(barcode_graph)))
    # print('max of barcode_graph {}'.format(statistics.mean((len(val) for val in barcode_graph.values()))))
    # for barcode, neighbors in barcode_graph.items():
    #     print(barcode, len(neighbors), neighbors)

    # log_file = open(args.log_file, 'a')
    # print('Building clusters', file=log_file)
    # log_file.close()
    # count = 0
    # barcode_clusters = [{i} for i in range(len(barcodes_1))]
    # a = time.time()
    # for barcode_set in lsh.values():
    #     count += 1
    #     if count % 100000 == 0:
    #         log_file = open(args.log_file, 'a')
    #         print('Count', count, file=log_file)
    #         log_file.close()
    #     union_set = set()
    #     for barcode in barcode_set:
    #         union_set.update(barcode_clusters[barcode])
    #     for barcode in union_set:
    #         barcode_clusters[barcode] = union_set
    #
    # b = time.time()
    # print("Time: ", b-a)
    #
    # id_dict = {}
    # for _set in barcode_clusters:
    #     id_dict[id(_set)] = _set
    # log_file = open(args.log_file, 'a')
    # print('Total number of clusters is:', len(id_dict), file=log_file)
    # print('len(cluster);', 'mean(d(v)) for v in cluster;', '[(v, edges(v)) for v in cluster]', file=log_file)
    # for _set in id_dict.values():
    #     if len(_set) > 0:
    #         print(len(_set),statistics.mean((len(barcode_graph[barcode]) for barcode in _set )),  [(barcode, barcode_graph[barcode]) for barcode in _set], sep='\t', file=log_file)
    # log_file.close()


if __name__ == '__main__':
    main()
