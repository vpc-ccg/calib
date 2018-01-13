import argparse

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
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    molecules_file = open(args.input_amplified_molecules)
    clusters_file = open(args.input_cluster_file)

    read_true_cluster = dict()
    true_cluster_reads = dict()
    predicted_clusters = list()

    for line in molecules_file.readlines():
        if line[0] != '>':
            continue
        line = line[1:].rstrip().split('_')
        read_id = line[0]
        cluster_id = line[3]
        if cluster_id in true_clusters:
            true_clusters[cluster_id].append(read_id)
        else:
            true_clusters[cluster_id] = [read_id]
        reads[read_id] = cluster_id

    for line in clusters_file.readlines():
        if line[0] != '#':
            predicted_clusters.append(list())
            continue
        line = line.split('\t')[2][1:].split('_')
        read_id = line[0]
        cluster_id = line[3]
        if cluster_id in true_clusters:
            true_clusters[cluster_id].append(read_id)
        else:
            true_clusters[cluster_id] = [read_id]







if __name__ == "__main__":
    main()
