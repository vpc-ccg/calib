import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculates Rand Index for a clusters file")
    parser.add_argument("-c",
                        "--input-clusters-file",
                        type=str,
                        required=True,
                        help="Input .clusters file to calculate Rand Index on")
    parser.add_argument("-m",
                        "--number-of-true-molecules",
                        type=int,
                        required=True,
                        help="The number of molecules generated in the simulated ground truth.")
    parser.add_argument("-r",
                        "--number-of-reads-per-true-molecule",
                        type=int,
                        required=True,
                        help="The number of simulted reads per molecule generated in the simulated ground truth.")
    args = parser.parse_args()
    return args


def score_cluster(cluster_counts):
    score = 0
    for i in range(len(cluster_counts)):
        score = score + cluster_counts[i]*(cluster_counts[i]-1)//2
        for j in range(i+1, len(cluster_counts)):
             score = score - cluster_counts[i]*cluster_counts[j]
    return score

def get_predicted_cluster(line):
    line = line.rstrip()
    line = line.split('\t')
    line = line[0]
    line = line.split('_')
    line = line[2]
    line = line.split(':')
    line = line[0]
    return int(line)

def main():
    args = parse_args()
    clusters_file = open(args.input_clusters_file)
    _molecules_total = args.number_of_true_molecules
    _reads_per_molecule = args.number_of_reads_per_true_molecule
    _reads_total = _molecules_total*_reads_per_molecule
    correct_clusters = 0
    incorrect_clusters = 0
    max_possible_TP = (_reads_total - _reads_per_molecule) * _reads_per_molecule * _molecules_total // 2
    line = clusters_file.readline()
    while(line != ''):
        line = clusters_file.readline()
        cluster_counts = dict()
        while (line != '' and line[0] != '#'):
            predicted_cluster = get_predicted_cluster(line)
            if predicted_cluster in cluster_counts:
                cluster_counts[predicted_cluster] = cluster_counts[predicted_cluster] + 1
            else:
                cluster_counts[predicted_cluster] = 1
            line = clusters_file.readline()
        cluster_counts = list(cluster_counts.values())
        if len(cluster_counts) == 1 and cluster_counts[0] == _reads_per_molecule:
            correct_clusters += 1
        else:
            incorrect_clusters += 1
        max_possible_TP = max_possible_TP + score_cluster(cluster_counts)
    print('Rand Index:', max_possible_TP/ (_reads_total*(_reads_total-1)//2 ))
    print('Correct clusters:', correct_clusters)
    print('Incorrect clutsers:', incorrect_clusters)


if __name__ == "__main__":
    main()
