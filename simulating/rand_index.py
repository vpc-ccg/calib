import sys
import argparse

from sklearn.metrics.cluster import adjusted_rand_score

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculates Rand Index for a clusters file")
    parser.add_argument("-p",
                        "--predicted-cluster-file",
                        type=str,
                        required=True,
                        help="Input predicted .cluster file")
    parser.add_argument("-t",
                        "--true-cluster-file",
                        type=str,
                        required=True,
                        help="Input true .cluster file")
    parser.add_argument("-o",
                        "--output-accuracy-results",
                        type=str,
                        required=False,
                        help="Output file where ARI score is outputted. Default: stdout")
    args = parser.parse_args()
    return args

def choose_2(n):
    return n*(n-1)//2

def main():
    args = parse_args()
    true_file = open(args.true_cluster_file)
    pred_file = open(args.predicted_cluster_file)

    n1_to_rid = dict()
    rid_to_tcid = dict()
    rid_to_pcid = dict()

    for line in true_file.readlines():
        line = line.rstrip().split('\t')
        cid = int(line[0])
        rid = int(line[2])
        n1  = line[3]

        n1_to_rid[n1] = rid
        rid_to_tcid[rid] = cid

    if (args.output_accuracy_results):
        results = open(args.output_accuracy_results, 'w+')
    else:
        results = sys.stdout

    for line in pred_file.readlines():
        line = line.rstrip().split('\t')
        cid = int(line[0])
        if line[3][0]=='@':
            line[3] = line[3][1:]
        if line[2] == "rid" and line[3] == "n1":
            print("ERR", line, file=results)
            exit()
        if line[2] != "rid" and line[3] == "n1":
            rid = int(line[2])
        if line[2] == "rid" and line[3] != "n1":
            rid = n1_to_rid[line[3]]
        if line[2] != "rid" and line[3] != "n1":
            rid = int(line[2])
            if rid != n1_to_rid[line[3]]:
                print("ERR", n1_to_rid[line[3]], line, file=results)
                exit()
        rid_to_pcid[rid] = cid
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

if __name__ == "__main__":
    main()
