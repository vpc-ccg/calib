import sys
import argparse
import subprocess

import plotly
from pymsa.score import Entropy

PARITITION_STRING = '{}: {}-{}'
OTHER_PARTITION_STRING = '{}: Other'

CID_WIDTH = 7
RID_WIDTH = 8

RNM_WIDTH = RID_WIDTH + 1 + 16
FLG_WIDTH = 3
CHR_WIDTH = 2
POS_WIDTH = 9
MPQ_WIDTH = 2
CIG_WIDTH = 16
FRG_WIDTH = 5
SEQ_WIDTH = 150
QUA_WIDTH = SEQ_WIDTH
TGS_WIDTH = 100

RNM_CLM = 0
FLG_CLM = 1
CHR_CLM = 2
POS_CLM = 3
MPQ_CLM = 4
CIG_CLM = 5
MCH_CLM = 6
MPS_CLM = 7
FRG_CLM = 8
SEQ_CLM = 9
QUA_CLM = 10
TGS_CLM = 11

RJST_CLM_TO_WIDTH = {
    RNM_CLM : RNM_WIDTH,
    FLG_CLM : FLG_WIDTH,
    CHR_CLM : CHR_WIDTH,
    POS_CLM : POS_WIDTH,
    MPQ_CLM : MPQ_WIDTH,
    MCH_CLM : CHR_WIDTH,
    MPS_CLM : POS_WIDTH,
    FRG_CLM : FRG_WIDTH
}
LJST_CLM_TO_WIDTH = {
    CIG_CLM : CIG_WIDTH,
    SEQ_CLM : SEQ_WIDTH,
    QUA_CLM : QUA_WIDTH
}
LEFT_POS_WIDTH = POS_WIDTH + CHR_WIDTH + 1
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
    parser.add_argument("-s",
                        "--sam-file",
                        type=str,
                        required=True,
                        help="Input sam file to extract mapping quality and mapping info from")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Output file where bipartite details are stored; Default: stdout")
    parser.add_argument("-w",
                        "--sankey-last-weight-index",
                        type=int,
                        required=False,
                        help="Python index of Sankey diagram's last layer's weight in the sorted list of (to be) observed weights; Default: min(|weights|-1, 4)")
    args = parser.parse_args()
    return args

def choose_2(n):
    return n*(n-1)//2

def main():
    args = parse_args()
    true_cluster_file = open(args.true_cluster_file)
    predicted_cluster_file = open(args.predicted_cluster_file)
    sam_file = open(args.sam_file)
    print('Good morning!')
    print('args.true_cluster_file', args.true_cluster_file, sep=':\t')
    print('args.predicted_cluster_file', args.predicted_cluster_file, sep=':\t')
    print('args.output', args.output, sep=':\t')
    print('args.sam_file', args.sam_file, sep=':\t')
    print('args.sankey_last_weight_index', args.sankey_last_weight_index, sep=':\t')

    reads = list()
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

    # contigency_matrix = dict()
    print('Reading SAM file')
    line = sam_file.readline()
    while len(line) > 0:
        if line[0] == '@':
            line = sam_file.readline()
            continue
        # rid = int(line.split('\t')[0].split('_')[0])
        line = line.rstrip().split('\t')
        other_line = sam_file.readline().rstrip().split('\t')
        for column, width in RJST_CLM_TO_WIDTH.items():
            line[column]       = line[column].rjust(width)
            other_line[column] = line[column].rjust(width)
        for column, width in LJST_CLM_TO_WIDTH.items():
            line[column]       = line[column].ljust(width)
            other_line[column] = line[column].ljust(width)
        line[TGS_CLM]       = ','.join(line[TGS_CLM:]      ).ljust(TGS_WIDTH)
        other_line[TGS_CLM] = ','.join(other_line[TGS_CLM:]).ljust(TGS_WIDTH)
        line = line[:TGS_CLM+1]
        other_line = other_line[:TGS_CLM+1]
        flag = int(line[1])
        if flag & 64:
            reads.append(line + other_line)
        else:
            reads.append(other_line + line)
        if len(reads) % 1e6 == 0:
            print('Read {} is read'.format(str(len(reads)).rjust(RID_WIDTH)))
            # print('\t'.join(reads[-1]))
        line = sam_file.readline()
        print('Read {} is read'.format(str(len(reads)).rjust(RID_WIDTH)))

    print('Reading true cluster file')
    tcid_counter = -1
    for line in true_cluster_file.readlines():
        if line[0] == '#':
            tcid_counter += 1
            tcid = tcid_counter
            if tcid % 1e5 == 1:
                print('Cluster {} is: {}'.format(str(tcid-1).rjust(CID_WIDTH), ','.join(str(rid).rjust(RID_WIDTH) for rid in tcid_to_rid_set[tcid - 1])))
            continue
        line = line.replace(' ','\t')
        line = line.rstrip().split('\t')
        rid = int(line[1])
        left_pos = str(line[0])
        reads[rid] = [str(tcid).rjust(CID_WIDTH), str(left_pos).rjust(LEFT_POS_WIDTH)] + reads[rid]
        if tcid in tcid_to_rid_set:
            tcid_to_rid_set[tcid].add(rid)
        else:
            tcid_to_rid_set[tcid] = {rid}
        rid_to_tcid[rid] = tcid
    print('Cluster {} is: {}'.format(str(tcid).rjust(CID_WIDTH), ','.join(str(rid).rjust(RID_WIDTH) for rid in tcid_to_rid_set[tcid])))

    print('Reading predicted cluster file')
    pcid_counter = -1
    for line in predicted_cluster_file.readlines():
        if line[0] == '#':
            pcid_counter += 1
            pcid = pcid_counter
            if pcid % 1e5 == 1:
                print('Cluster {} is: {}'.format(str(pcid-1).rjust(CID_WIDTH), ','.join(str(rid).rjust(RID_WIDTH) for rid in pcid_to_rid_set[pcid - 1])))
            continue
        line = line.replace(' ', '\t')
        line = line.rstrip().split('\t')
        rid = int(line[1])
        tcid = rid_to_tcid[rid]
        rid_to_pcid[rid] = pcid

        minimizers_1 = '-'.join(line[3].split('-')[1:])
        minimizers_2 = '-'.join(line[6].split('-')[1:])
        reads[rid] = [str(rid).rjust(RID_WIDTH), str(pcid).rjust(CID_WIDTH), minimizers_1, minimizers_2] + reads[rid]

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
    print('Cluster {} is: {}'.format(str(pcid).rjust(CID_WIDTH), ','.join(str(rid).rjust(RID_WIDTH) for rid in pcid_to_rid_set[pcid])))

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

    weights = set()
    pc_size = len(pcidx_to_pcid)
    tc_size = len(tcidx_to_tcid)

    adjacency = [None]*node_count

    print('Building bipartite graph')
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
            neighbors[pcidx] = len(tcid_to_rid_set[tcidx_to_tcid[tcidx]].intersection(pcid_to_rid_set[pcidx_to_pcid[pcidx]]))
            weights.add(weight)
        adjacency[tcidx] = neighbors

    sankey_output = args.output

    weights = sorted(weights)
    history = [(0,0)]*node_count
    prev_weight = -1

    if args.sankey_last_weight_index and args.sankey_last_weight_index >= -len(weights) and args.sankey_last_weight_index < len(weights):
        sankey_last_weight = weights[args.sankey_last_weight_index]
    else:
        if 4 < len(weights):
            sankey_last_weight = weights[4]
        else:
            sankey_last_weight = weights[-1]


    sankey_labels = list()
    sankey_colors = list()
    # sankey_labels.append(PARITITION_STRING.format(weights[0],  1 ,  1 ))
    # sankey_colors.append('black')
    # sankey_labels.append(PARITITION_STRING.format(weights[0],  1 , 'X'))
    # sankey_colors.append('#bdd7e7')
    # sankey_labels.append(PARITITION_STRING.format(weights[0],  1 ,  2 ))
    # sankey_colors.append('#6baed6')
    # sankey_labels.append(PARITITION_STRING.format(weights[0],  1 ,  0 ))
    # sankey_colors.append('#2171b5')
    # sankey_labels.append(PARITITION_STRING.format(weights[0], 'X',  1 ))
    # sankey_colors.append('#bae4b3')
    # sankey_labels.append(PARITITION_STRING.format(weights[0],  2 ,  1 ))
    # sankey_colors.append('#74c476')
    # sankey_labels.append(PARITITION_STRING.format(weights[0],  0 ,  1 ))
    # sankey_colors.append('#238b45')
    # sankey_labels.append(OTHER_PARTITION_STRING.format(weights[0]))
    # sankey_colors.append('white')
    for weight in weights:
        if weight > sankey_last_weight:
            break
        sankey_labels.append(PARITITION_STRING.format(weight,  1 ,  1 ))
        sankey_colors.append('black')
        sankey_labels.append(PARITITION_STRING.format(weight,  1 , 'X'))
        sankey_colors.append('#bdd7e7')
        sankey_labels.append(PARITITION_STRING.format(weight,  1 ,  0 ))
        sankey_colors.append('#2171b5')
        sankey_labels.append(PARITITION_STRING.format(weight, 'X',  1 ))
        sankey_colors.append('#bae4b3')
        sankey_labels.append(PARITITION_STRING.format(weight,  0 ,  1 ))
        sankey_colors.append('#238b45')
        sankey_labels.append(OTHER_PARTITION_STRING.format(weight))
        sankey_colors.append('white')
    sankey_label_to_id = dict()
    for idx, label in enumerate(sankey_labels):
        sankey_label_to_id[label] = idx
    sankey_read_link_dict = dict()
    sankey_clus_link_dict = dict()

    for current_weight in weights[:sankey_last_weight]:
        print('Building sankey diagram column with weight {}'.format(current_weight))
        visited = [False]*node_count
        shared_cc = 0
        pc_cc = 0
        tc_cc = 0
        cc_composition = dict()
        cc_examples = dict()
        for idx in range(node_count):
            if visited[idx]:
                continue
            cc_list = list()
            visited[idx] = True
            cc_list.append(idx)
            opened = [idx]
            while len(opened) > 0:
                current_node = opened.pop()
                for neighbor in adjacency[current_node]:
                    weight = adjacency[current_node][neighbor]
                    if weight >= current_weight and visited[neighbor] != True:
                        visited[neighbor] = True
                        cc_list.append(neighbor)
                        opened.append(neighbor)
            cc_rid_list = list()
            cc_pc_count = 0
            cc_tc_count = 0
            for current_node in cc_list:
                if current_node < pc_size:
                    cc_rid_list.extend(pcid_to_rid_set[pcidx_to_pcid[current_node]])
                    cc_pc_count += 1
                else:
                    cc_rid_list.extend(tcid_to_rid_set[tcidx_to_tcid[current_node]])
                    cc_tc_count += 1
            key = (cc_pc_count, cc_tc_count, len(cc_rid_list))
            if key in cc_composition:
                cc_composition[key] += 1
            else:
                cc_composition[key] = 1
                cc_examples[key] = '\n'.join(('\t'.join(reads[rid]) for rid in set(cc_rid_list)))
            for current_node in cc_list:
                if current_weight > sankey_last_weight:
                    break
                last_cc_pc_count = history[current_node][0]
                last_cc_tc_count = history[current_node][1]
                history[current_node] = (cc_pc_count, cc_tc_count)
                if prev_weight == -1:
                    continue
                source = PARITITION_STRING.format(prev_weight, last_cc_pc_count, last_cc_tc_count)
                if source not in sankey_label_to_id:
                    if   last_cc_pc_count == 1 and last_cc_tc_count > 1:
                        source = PARITITION_STRING.format(prev_weight, 1, 'X')
                    elif last_cc_pc_count > 1 and last_cc_tc_count == 1:
                        source = PARITITION_STRING.format(prev_weight, 'X', 1)
                    else:
                        source = OTHER_PARTITION_STRING.format(prev_weight)
                source_id = sankey_label_to_id[source]

                target = PARITITION_STRING.format(current_weight, cc_pc_count, cc_tc_count)
                if target not in sankey_label_to_id:
                    if   cc_pc_count == 1 and cc_tc_count > 1:
                        target = PARITITION_STRING.format(current_weight, 1, 'X')
                    elif cc_pc_count > 1 and cc_tc_count == 1:
                        target = PARITITION_STRING.format(current_weight, 'X', 1)
                    else:
                        target = OTHER_PARTITION_STRING.format(current_weight)
                target_id = sankey_label_to_id[target]

                key = (source_id, target_id)

                if current_node < pc_size:
                    node_rid_count = len(pcid_to_rid_set[pcidx_to_pcid[current_node]])
                else:
                    node_rid_count = len(tcid_to_rid_set[tcidx_to_tcid[current_node]])

                if key in sankey_clus_link_dict:
                    sankey_clus_link_dict[key] += 1
                else:
                    sankey_clus_link_dict[key] = 1

                if key in sankey_read_link_dict:
                    sankey_read_link_dict[key] += node_rid_count
                else:
                    sankey_read_link_dict[key] = node_rid_count
        output = open(args.output+'.w_'+str(current_weight)+'.log', 'w+')
        for key in sorted(cc_composition):
            cc_pc_count = int(key[0])
            cc_tc_count = int(key[1])
            cc_rid_list_len = int(key[2])
            cc_count = int(cc_composition[key])
            print('#',int(current_weight), cc_pc_count, cc_tc_count, cc_rid_list_len, cc_count, sep='\t', file=output)
            print(cc_examples[key], file=output)
        output.close()
        prev_weight = current_weight

    data = dict(
        type='sankey',
        node = dict(
          pad = 15,
          thickness = 20,
          line = dict(
            color = "black",
            width = 0.5
          ),
          label = sankey_labels,
          color = sankey_colors
        ))

    layout =  dict(
        font = dict(
          size = 10
        )
    )

    sources = list()
    targets = list()
    values = list()
    for key in sankey_read_link_dict:
        sources.append(key[0])
        targets.append(key[1])
        values.append(sankey_read_link_dict[key])
    data['link'] = dict(source = sources,  target = targets, value = values)
    layout['title'] = "Sankey diagram of reads"
    fig = dict(data=[data], layout=layout)
    plotly.offline.plot(fig, filename=sankey_output+'.sankey_reads.html', auto_open=False)

    sources = list()
    targets = list()
    values = list()
    for key in sankey_clus_link_dict:
        sources.append(key[0])
        targets.append(key[1])
        values.append(sankey_clus_link_dict[key])
    data['link'] = dict(source = sources,  target = targets, value = values)
    layout['title'] = "Sankey diagram of clusters"
    fig = dict(data=[data], layout=layout)
    plotly.offline.plot(fig, filename=sankey_output+'.sankey_clusters.html', auto_open=False)

if __name__ == "__main__":
    main()
