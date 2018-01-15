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
    clustering_reads = list()

    for line in molecules_file.readlines():
        if line[0] != '>':
            continue
        line = line[1:].rstrip().split('_')
        read_id = line[0]
        cluster_id = line[3]
        if cluster_id in true_cluster_reads:
            true_cluster_reads[cluster_id].append(read_id)
        else:
            true_cluster_reads[cluster_id] = [read_id]
        read_true_cluster[read_id] = cluster_id


    contigency_matrix = dict()
    predicted_cluster_count = 0

    true_cluster_sum = 0
    predicted_cluster_sum = 0
    for line in clusters_file.readlines():
        if line[0] != '#':
            predicted_cluster_count += 1
            continue
        line = line.split('\t')[2][1:].split('_')
        # read_id = line[0]
        true_cluster_id = line[3]
        key = (true_cluster_id, predicted_cluster_count-1)
        contigency_matrix[key] = contigency_matrix.get(key, 0) + 1

    rand_index = 0
    for value in contigency_matrix.values():
        rand_index += value*(value-1)/2


    true_cluster_sum = [0 for _ in range(len(true_cluster_reads))]
    predicted_cluster_sum = [0 for_ in range(predicted_cluster_count)]







if __name__ == "__main__":
    main()
