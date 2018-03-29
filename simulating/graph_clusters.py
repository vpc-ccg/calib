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
        reads[rid] = str(line[0])
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
        tcid = rid_to_tcid[rid]

        rid_to_pcid[rid] = pcid_counter
        reads[rid] = line[3]+'\t'+line[6]+'\t'+str(rid)+'\t'+reads[rid]+'\t'+str(rid_to_pcid[rid])+'\t'+str(rid_to_tcid[rid])
        # print(reads[rid])
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

    # if (args.output_accuracy_results):
    #     results = open(args.output_accuracy_results+'.ari', 'w+')
    # else:
    #     results = sys.stdout

    # Y = []
    # X = []
    # lonely_read_cluster_counter = -1
    # for rid in rid_to_tcid:
    #     Y.append(rid_to_tcid[rid])
    #     if rid in rid_to_pcid:
    #         X.append(rid_to_pcid[rid])
    #     else:
    #         X.append(lonely_read_cluster_counter)
    #         lonely_read_cluster_counter -= 1
    pcid_to_pcidx = dict()
    tcid_to_tcidx = dict()
    pcidx_to_pcid = dict()
    tcidx_to_tcid = dict()
    node_count = 0

    pcid_singltons = set()
    tcid_singltons = set()
    for rid in rid_to_tcid:
        pcid = rid_to_pcid[rid]
        tcid = rid_to_tcid[rid]
        if len(tcid_to_rid_set[tcid]) == 1 and len(pcid_to_rid_set[pcid]) == 1:
            tcid_singltons.add(tcid)
            pcid_singltons.add(pcid)

    for pcid in pcid_to_tcid_set:
        # if pcid in pcid_singltons:
        #     pass
        pcidx_to_pcid[node_count] = pcid
        pcid_to_pcidx[pcid] = node_count
        node_count += 1
    for tcid in tcid_to_pcid_set:
        # if tcid in tcid_singltons:
        #     pass
        tcidx_to_tcid[node_count] = tcid
        tcid_to_tcidx[tcid] = node_count
        node_count += 1

    weights = set()
    pc_size = len(pcidx_to_pcid)
    tc_size = len(tcidx_to_tcid)

    adjacency = [None]*node_count

    for pcidx in range(pc_size):
        neighbors = {tcid_to_tcidx[tcid] : 0 for tcid in pcid_to_tcid_set[pcidx_to_pcid[pcidx]]}
        for tcidx in neighbors:
            weight = len(pcid_to_rid_set[pcidx_to_pcid[pcidx]].intersection(tcid_to_rid_set[tcidx_to_tcid[tcidx]]))
            neighbors[tcidx] = weight
            weights.add(weight)
        adjacency[pcidx] = neighbors
    for tcidx in range(pc_size, node_count):
        neighbors = {pcid_to_pcidx[pcid] : 0 for pcid in tcid_to_pcid_set[tcidx_to_tcid[tcidx]]}
        for pcidx in neighbors:
            # weight =
            neighbors[pcidx] = len(tcid_to_rid_set[tcidx_to_tcid[tcidx]].intersection(pcid_to_rid_set[pcidx_to_pcid[pcidx]]))
            weights.add(weight)
        adjacency[tcidx] = neighbors

    if (args.output_accuracy_results):
        output = open(args.output_accuracy_results, 'w+')
    else:
        output = sys.stdout

    for current_weight in sorted((weights)):
        visited = [False]*node_count
        shared_cc = 0
        pc_cc = 0
        tc_cc = 0
        cc_composition = dict()
        cc_examples = dict()
        for idx in range(node_count):
            # print (idx, visited[idx])
            if visited[idx]:
                continue
            cc_size = 0
            cc_rid_list = list()
            cc_pc_count = 0
            cc_tc_count = 0
            visited[idx] = True
            opened = [idx]
            # print(opened)
            while len(opened) > 0:
                # print (idx, opened)
                current_node = opened.pop()
                if current_node < pc_size:
                    cc_rid_list.extend(list(pcid_to_rid_set[pcidx_to_pcid[current_node]]))
                    cc_pc_count += 1
                else:
                    cc_rid_list.extend(list(tcid_to_rid_set[tcidx_to_tcid[current_node]]))
                    cc_tc_count += 1
                # print(opened)
                cc_size += 1

                # print (adjacency[current_node], type(adjacency[current_node]))
                for neighbor in adjacency[current_node]:
                    weight = adjacency[current_node][neighbor]
                    # print ('\t-->', neighbor, weight)
                    if weight >= current_weight and visited[neighbor] != True:
                        visited[neighbor] = True
                        opened.append(neighbor)
            if cc_size ==  1:
                if idx < pc_size:
                    pc_cc += 1
                else:
                    tc_cc += 1
            else:
                shared_cc += 1
            rid_count = 0
            cc_rid_list.sort()
            cc_rid_list_unique = list()
            # print(cc_rid_list)
            for index, rid in enumerate(cc_rid_list):
                if index == len(cc_rid_list) - 1 or rid != cc_rid_list[index+1]:
                    cc_rid_list_unique.append(rid)
            # print(rid_count)
            key = (cc_pc_count, cc_tc_count, len(cc_rid_list_unique))
            if key in cc_composition:
                cc_composition[key] += 1
            else:
                cc_composition[key] = 1
                cc_examples[key] = '\n'.join((reads[rid] for rid in cc_rid_list_unique))


        for key in sorted(cc_composition):
            print('#',int(current_weight), int(key[0]), int(key[1]), int(key[2]), int(cc_composition[key]), sep='\t', file=output)
            print(cc_examples[key], file=output)

if __name__ == "__main__":
    main()
