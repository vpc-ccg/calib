import sys
import argparse

from sklearn.metrics.cluster import adjusted_rand_score

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
    print(adjusted_rand_score(X, Y), file=results)

    if (args.output_accuracy_results):
        results = open(args.output_accuracy_results+'.merged', 'w+')
    else:
        results = sys.stdout

    # print('Clusters merged in prediction:\n===', file=results)
    for pcid in pcid_to_tcid_set:
        if len(pcid_to_tcid_set[pcid]) > 1:
            true_clusters = [(x, len(tcid_to_rid_set[x])) for x in pcid_to_tcid_set[pcid]]
            # print(pcid,'\n', true_clusters)
            print('#\t{} have reads merged into {}'.format(true_clusters, tuple([pcid, len(pcid_to_rid_set[pcid])]) ), file=results)
            for tcid in pcid_to_tcid_set[pcid]:
                for rid in tcid_to_rid_set[tcid].intersection(pcid_to_rid_set[pcid]):
                    print('{}\t{}\t{}\t{}'.format(reads[rid], tcid, rid_to_pcid[rid], rid), file=results)
                # print(file=results)

    if (args.output_accuracy_results):
        results = open(args.output_accuracy_results+'.split', 'w+')
    else:
        results = sys.stdout

    # print('Clusters split in prediction:\n===', file=results)
    for tcid in tcid_to_pcid_set:
        if len(tcid_to_pcid_set[tcid]) > 1:
            predicted_clusters = [(x, len(pcid_to_rid_set[x])) for x in tcid_to_pcid_set[tcid]]
            print('#\t{} have reads split into {}'.format(tuple([tcid, len(tcid_to_rid_set[tcid])]), predicted_clusters), file=results)
            for pcid in tcid_to_pcid_set[tcid]:
                for rid in pcid_to_rid_set[pcid].intersection(tcid_to_rid_set[tcid]):
                    print('{}\t{}\t{}\t{}'.format(reads[rid], tcid, rid_to_pcid[rid], rid), file=results)


    # for tcid in tcid_to_pcid_set:
    #     if len(tcid_to_pcid_set[tcid]) > 1:
    #         print('|', tcid, '| =', len(tcid_to_rid_set[tcid]), file=results)
    #         for pcid in tcid_to_pcid_set[tcid]:
    #             print('\t|', tcid, '∩', pcid,'| =', len(pcid_to_rid_set[pcid].intersection(tcid_to_rid_set[tcid])), file=results)
    #




if __name__ == "__main__":
    main()
