import sys
import argparse

from sklearn.metrics.cluster import adjusted_rand_score

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculates Rand Index for a clusters file")
    parser.add_argument("-c",
                        "--input-cluster-file",
                        type=str,
                        required=True,
                        help="Input .cluster file to calculate Rand Index on")
    parser.add_argument("-m",
                        "--input-amplified-molecules",
                        type=str,
                        required=True,
                        help="Input .fa file of .")
    parser.add_argument("-o",
                        "--output-accuracy-results",
                        type=str,
                        required=False,
                        help="Output file where accuracy results and contents of any clutsters with clustering mistakes will be printed. Default: stdout")
    args = parser.parse_args()
    return args

def choose_2(n):
    return n*(n-1)//2

def main():
    args = parse_args()
    molecules_file = open(args.input_amplified_molecules)
    clusters_file = open(args.input_cluster_file)

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

    for line in molecules_file.readlines():
        if line[0] != '>':
            continue
        line = line[1:].rstrip().split('_')
        rid = int(line[0])
        tcid = int(line[3].split(':')[0])

        if tcid in tcid_to_rid_set:
            tcid_to_rid_set[tcid].add(rid)
        else:
            tcid_to_rid_set[tcid] = {rid}
        rid_to_tcid[rid] = tcid

    pcid_counter = -1
    for line in clusters_file.readlines():
        if line[0] == '#':
            pcid_counter += 1
            continue

        line = line.split('\t')[2][1:].split('_')
        rid = int(line[0])
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
        results = open(args.output_accuracy_results, 'w+')
    else:
        results = sys.stdout

    Y = []
    X = []
    for rid in rid_to_tcid:
        Y.append(rid_to_tcid[rid])
        X.append(rid_to_pcid[rid])
    print(adjusted_rand_score(X, Y), file=results)

    for pcid in pcid_to_tcid_set:
        if len(pcid_to_tcid_set[pcid]) > 1:
            print(pcid, file=results)
            for tcid in pcid_to_tcid_set[pcid]:
                print('\t{}'.format(tcid), end='', file=results)
                for rid in tcid_to_rid_set[tcid].intersection(pcid_to_rid_set[pcid]):
                    print('\t{}'.format(rid), end='', file=results)
                print(file=results)

    print('=', file=results)
    for tcid in tcid_to_pcid_set:
        if len(tcid_to_pcid_set[tcid]) > 1:
            print(tcid, file=results)
            for pcid in tcid_to_pcid_set[tcid]:
                print('\t{}'.format(pcid), end='', file=results)
                for rid in pcid_to_rid_set[pcid].intersection(tcid_to_rid_set[tcid]):
                    print('\t{}'.format(rid), end='', file=results)
                print(file=results)


    # for tcid in tcid_to_pcid_set:
    #     if len(tcid_to_pcid_set[tcid]) > 1:
    #         print('|', tcid, '| =', len(tcid_to_rid_set[tcid]), file=results)
    #         for pcid in tcid_to_pcid_set[tcid]:
    #             print('\t|', tcid, '∩', pcid,'| =', len(pcid_to_rid_set[pcid].intersection(tcid_to_rid_set[tcid])), file=results)
    #




if __name__ == "__main__":
    main()
