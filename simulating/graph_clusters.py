import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculates Rand Index for a clusters file")
    parser.add_argument("-x",
                        "--predicted-cluster-file",
                        type=str,
                        required=True,
                        help="Input .cluster file to calculate Rand Index on")
    parser.add_argument("-y",
                        "--true-cluster-file",
                        type=str,
                        required=True,
                        help="Input .cluster file to calculate Rand Index on")
    parser.add_argument("-o",
                        "--output-accuracy-results",
                        type=str,
                        required=False,
                        help="Output file where accuracy results and contents of any clusters with clustering discordances will be printed. Default: stdout")
    args = parser.parse_args()
    return args

def choose_2(n):
    return n*(n-1)//2

def main():
    args = parse_args()
    predicted_cluster_file = open(args.predicted_cluster_file)
    true_cluster_file = open(args.true_cluster_file)

    reads = dict()
    # rid_to_tcid --> read_id_to_true_cluster_id
    rid_to_tcid = dict()
    # rid_to_pcid --> read_id_to_predicted_cluster_id
    rid_to_pcid = dict()
    # tcid_to_rid_set --> true_cluster_id_to_read_ids_set
    tcid_to_rid_set = dict()
    # pcid_to_rid_set --> predicted_cluster_id_to_read_ids_set
    pcid_to_rid_set = dict()
    # pcid_to_tcid_set --> predicted_cluster_id_to_true_cluster_ids_set
    pcid_to_tcid_set = dict()
    # tcid_to_pcid_set --> true_cluster_id_to_predicted_cluster_ids_set
    tcid_to_pcid_set = dict()

    contigency_matrix = dict()

    tcid_counter = -1
    for line in true_cluster_file.readlines():
        if line[0] == '#':
            tcid_counter += 1
            continue
        line = line.replace(' ','\t')
        line = line.rstrip().split('\t')
        tcid = tcid_counter
        rid = int(line[1])

        if tcid in tcid_to_rid_set:
            tcid_to_rid_set[tcid].add(rid)
        else:
            tcid_to_rid_set[tcid] = {rid}
        rid_to_tcid[rid] = tcid

    pcid_counter = -1
    for line in predicted_cluster_file.readlines():
        if line[0] == '#':
            pcid_counter += 1
            continue
        line = line.replace(' ', '\t')
        line = line.rstrip().split('\t')
        rid = int(line[1])
        reads[rid] = line[3]+'\t'+line[6]
        tcid = rid_to_tcid[rid]

        rid_to_pcid[rid] = pcid_counter
        if pcid_counter in pcid_to_rid_set:
            pcid_to_rid_set[pcid_counter].add(rid)
        else:
            pcid_to_rid_set[pcid_counter] = {rid}

        if pcid_counter in pcid_to_tcid_set:
            pcid_to_tcid_set[pcid_counter].add(tcid)
        else:
            pcid_to_tcid_set[pcid_counter] = {tcid}

        if tcid in tcid_to_pcid_set:
            tcid_to_pcid_set[tcid].add(pcid_counter)
        else:
            tcid_to_pcid_set[tcid] = {pcid_counter}

        key = tuple([tcid, pcid_counter])
        contigency_matrix[key] = contigency_matrix.get(key, 0) + 1

    if (args.output_accuracy_results):
        results = open(args.output_accuracy_results+'.ari', 'w+')
    else:
        results = sys.stdout

    Y = []
    X = []
    lonely_read_cluster_counter = -1
    for rid in rid_to_tcid:
        Y.append(rid_to_tcid[rid])
        if rid in rid_to_pcid:
            X.append(rid_to_pcid[rid])
        else:
            X.append(lonely_read_cluster_counter)
            lonely_read_cluster_counter -= 1
    pcid_to_pcidx = dict()
    tcid_to_tcidx = dict()
    pcidx_to_pcid = dict()
    tcidx_to_tcid = dict()
    node_count = 0
    for pcid in pcid_to_tcid_set:
        pcidx_to_pcid[node_count] = pcid
        pcid_to_pcidx[pcid] = node_count
        node_count += 1
    for tcid in tcid_to_pcid_set:
        tcidx_to_tcid[node_count] = tcid
        tcid_to_tcidx[tcid] = node_count
        node_count += 1

    max_weight = 0
    pc_size = len(pcid_to_tcid_set)
    tc_size = len(tcid_to_pcid_set)

    adjacency = [None]*node_count

    for pcidx in range(pc_size):
        neighbors = {tcid_to_tcidx[tcid] : 0 for tcid in pcid_to_tcid_set[pcidx_to_pcid[pcidx]]}
        for tcidx in neighbors:
            weight = len(pcid_to_rid_set[pcidx_to_pcid[pcidx]].intersection(tcid_to_rid_set[tcidx_to_tcid[tcidx]]))
            neighbors[tcidx] = weight
            if weight > max_weight:
                max_weight = weight
        adjacency[pcidx] = neighbors
    for tcidx in range(pc_size, node_count):
        neighbors = {pcid_to_pcidx[pcid] : 0 for pcid in tcid_to_pcid_set[tcidx_to_tcid[tcidx]]}
        for pcidx in neighbors:
            # weight =
            neighbors[pcidx] = len(tcid_to_rid_set[tcidx_to_tcid[tcidx]].intersection(pcid_to_rid_set[pcidx_to_pcid[pcidx]]))
            if neighbors[pcidx] > max_weight:
                max_weight = neighbors[pcidx]
        adjacency[tcidx] = neighbors

    for current_weight in range(max_weight+2):
        visited = [False]*node_count
        shared_cc = 0
        pc_cc = 0
        tc_cc = 0
        for idx in range(node_count):
            # print (idx, visited[idx])
            if visited[idx]:
                continue
            cc_size = 1
            visited[idx] = True
            opened = [idx]
            # print(opened)
            while len(opened) > 0:
                # print (opened)
                current_node = opened.pop()
                # print (adjacency[current_node], type(adjacency[current_node]))
                for neighbor in adjacency[current_node]:
                    weight = adjacency[current_node][neighbor]
                    # print ('\t-->', neighbor, weight)
                    if weight >= current_weight and visited[neighbor] != True:
                        visited[neighbor] = True
                        opened.append(neighbor)
                        cc_size += 1
            if cc_size ==  1:
                if idx < pc_size:
                    pc_cc += 1
                else:
                    tc_cc += 1
            else:
                shared_cc += 1
        print(current_weight, pc_cc, tc_cc, shared_cc, sep='\t')

if __name__ == "__main__":
    main()
