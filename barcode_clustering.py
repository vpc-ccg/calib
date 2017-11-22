import sys
import time
import statistics
from itertools import combinations
_complements = {'A': 'T',
                'T': 'A',
                'G': 'C',
                'C': 'G'}


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
    nodes_to_traverse = graph[node] - cluster
    cluster.update(graph[node])
    cluster.add(node)
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


def main():
    _barcode_length = 10
    _error_tolerance = 2
    log_file = open(sys.argv[3], 'w+')
    print('Hmmmm. Good morning?!', file=log_file)
    log_file.close()
    barcode_lines_1 = get_barcodes(sys.argv[1], _barcode_length)
    barcode_lines_2 = get_barcodes(sys.argv[2], _barcode_length)
    if not len(barcode_lines_1) == len(barcode_lines_2):
        log_file = open(sys.argv[3], 'a')
        print('You messed up; the read files are not of the same length', file=log_file)
        log_file.close()
        sys.exit(42)
    log_file = open(sys.argv[3], 'a')
    print('Let\'s start this shit!', file=log_file)
    barcode_pairs_to_lines = get_barcode_pair_to_line_mapping(barcode_lines_1, barcode_lines_2)

    print('Total number of unique barcode pairs:', len(barcode_pairs_to_lines), file=log_file)

    # for p in barcode_pairs_to_lines:
    #     print(p, len(barcode_pairs_to_lines[p]), barcode_pairs_to_lines[p])

    barcodes_1 = ['']*len(barcode_pairs_to_lines)
    barcodes_2 = ['']*len(barcode_pairs_to_lines)
    barcodes_rev_compl_1 = ['']*len(barcode_pairs_to_lines)
    barcodes_rev_compl_2 = ['']*len(barcode_pairs_to_lines)

    for index, barcode_pair in enumerate(barcode_pairs_to_lines.keys()):
        barcodes_1[index] = barcode_pair[0]
        barcodes_2[index] = barcode_pair[1]
        barcodes_rev_compl_1[index] = reverse_complement(barcode_pair[0])
        barcodes_rev_compl_2[index] = reverse_complement(barcode_pair[1])

    log_file.close()
    lsh = {}
    for template, template_id in template_generator(_barcode_length, _error_tolerance):
        # log_file = open(sys.argv[3], 'a')
        # print("Template {} with ID {}".format(template('12345678'), template_id), file=log_file)
        # log_file.close()
        for barcode_num in range(len(barcodes_1)):
            barcode_1 = template(barcodes_1[barcode_num])
            barcode_1_rev = template(barcodes_rev_compl_1[barcode_num])
            barcode_2 = template(barcodes_2[barcode_num])
            barcode_2_rev = template(barcodes_rev_compl_2[barcode_num])

            new_key = barcode_1 + barcode_2 + str(template_id)
            lsh[new_key] = lsh.get(new_key, {barcode_num}).union({barcode_num})
            new_key = barcode_2_rev + barcode_1_rev + str(template_id)
            lsh[new_key] = lsh.get(new_key, {barcode_num}).union({barcode_num})


    barcode_graph = dict()
    log_file = open(sys.argv[3], 'a')
    print("Building barcode graph", file=log_file)
    log_file.close()
    barcode_pair_length = _barcode_length*2 - _error_tolerance*2
    count = 0
    for val in lsh.values():
        count += 1
        if count % 100000 == 0:
            log_file = open(sys.argv[3], 'a')
            print('Count', count, file=log_file)
            log_file.close()
        for node in val:
            adjacent_nodes = val#.difference({node})
            if node in barcode_graph:
                barcode_graph[node].update(adjacent_nodes)
            else:
                barcode_graph[node] = adjacent_nodes
    a = time.time()
    clusters = generate_clusters_by_bfs(barcode_graph)
    b = time.time()
    print(b - a)
    print(len(clusters))
    print([len(cluster) for cluster in clusters])

    # print('len of barcode_graph {}'.format(len(barcode_graph)))
    # print('max of barcode_graph {}'.format(statistics.mean((len(val) for val in barcode_graph.values()))))
    # for barcode, neighbors in barcode_graph.items():
    #     print(barcode, len(neighbors), neighbors)

    log_file = open(sys.argv[3], 'a')
    print('Building clusters', file=log_file)
    log_file.close()
    count = 0
    barcode_clusters = [{i} for i in range(len(barcodes_1))]
    a = time.time()
    for barcode_set in lsh.values():
        count += 1
        if count % 100000 == 0:
            log_file = open(sys.argv[3], 'a')
            print('Count', count, file=log_file)
            log_file.close()
        union_set = set()
        for barcode in barcode_set:
            union_set.update(barcode_clusters[barcode])
        for barcode in union_set:
            barcode_clusters[barcode] = union_set

    b = time.time()
    print("Time: ", b-a)

    id_dict = {}
    for _set in barcode_clusters:
        id_dict[id(_set)] = _set
    log_file = open(sys.argv[3], 'a')
    print('Total number of clusters is:', len(id_dict), file=log_file)
    print('len(cluster);', 'mean(d(v)) for v in cluster;', '[(v, edges(v)) for v in cluster]', file=log_file)
    for _set in id_dict.values():
        if len(_set) > 0:
            print(len(_set),statistics.mean((len(barcode_graph[barcode]) for barcode in _set )),  [(barcode, barcode_graph[barcode]) for barcode in _set], sep='\t', file=log_file)
    log_file.close()


if __name__ == '__main__':
    main()
