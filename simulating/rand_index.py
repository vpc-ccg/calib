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
    args = parser.parse_args()
    return args

def choose_2(n):
    return n*(n-1)//2

def main():
    args = parse_args()
    molecules_file = open(args.input_amplified_molecules)
    clusters_file = open(args.input_cluster_file)

    read_id_to_true_cluster_id = dict()
    read_id_to_predicted_cluster_id = dict()

    true_cluster_id_to_read_id_list = dict()
    predicted_cluster_id_to_read_id_list = list()

    for line in molecules_file.readlines():
        if line[0] != '>':
            continue
        line = line[1:].rstrip().split('_')
        read_id = line[0]
        cluster_id = line[3].split(':')[0]
        # print(read_id, cluster_id)
        if cluster_id in true_cluster_id_to_read_id_list:
            true_cluster_id_to_read_id_list[cluster_id].append(read_id)
        else:
            true_cluster_id_to_read_id_list[cluster_id] = [read_id]
        read_id_to_true_cluster_id[read_id] = cluster_id
    read_count = len(read_id_to_true_cluster_id)
    # print('Read count is {}'.format(read_count))

    contigency_matrix = dict()
    predicted_cluster_count = 0

    true_cluster_sum = 0
    predicted_cluster_sum = 0
    for line in clusters_file.readlines():
        if line[0] == '#':
            predicted_cluster_id_to_read_id_list.append(list())
            predicted_cluster_count += 1
            continue
        line = line.split('\t')[2][1:].split('_')
        read_id = line[0]
        read_id_to_predicted_cluster_id[read_id]=predicted_cluster_count

        predicted_cluster_id_to_read_id_list[-1].append(read_id)
        true_cluster_id = int(line[3].split(':')[0])
        # print(predicted_cluster_count, true_cluster_id)

        key = tuple([true_cluster_id, predicted_cluster_count-1])
        contigency_matrix[key] = contigency_matrix.get(key, 0) + 1
        # print(key, contigency_matrix[key])
    # print (len(contigency_matrix))
    # for k in contigency_matrix:
    #     print(k, contigency_matrix[k])
    A_sum = sum((choose_2(len(a)) for a in true_cluster_id_to_read_id_list))
    B_sum = sum((choose_2(len(b)) for b in predicted_cluster_id_to_read_id_list))
    B = [len(b) for b in predicted_cluster_id_to_read_id_list]
    I = [n for n in contigency_matrix.values()]

    # index = sum((choose_2(n) for n in contigency_matrix.values()))
    # expected_index = A_sum*B_sum/choose_2(read_count)
    # max_index = (A_sum + B_sum)/2
    # ARI = (index - expected_index)/(max_index - expected_index)
    # print('A_sum {}\nB_sum {}\nindex {}\nexpected_index {}\nmax_index {}\nARI {}'.format(A_sum, B_sum, index, expected_index, max_index, ARI))

    Y = []
    X = []
    for read in read_id_to_true_cluster_id:
        Y.append(read_id_to_true_cluster_id[read])
        X.append(read_id_to_predicted_cluster_id[read])
    print(adjusted_rand_score(X, Y))

    # f = open('test', 'w+')
    # for k in contigency_matrix:
    #     print('{}\t{}\t{}'.format(k[0],k[1], contigency_matrix[k]) ,file=f)
    # f = open('sett1', 'w+')
    # for k in sorted(B):
    #     print(k ,file=f)
    # f = open('sett2', 'w+')
    # for k in sorted(I):
    #     print(k ,file=f)






if __name__ == "__main__":
    main()
