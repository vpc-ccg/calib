import argparse

import networkx as nx

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculates Rand Index for a clusters file")
    parser.add_argument("-i",
                        "--input-debug-cluster",
                        type=str,
                        required=True,
                        help='input debug cluster file')
    parser.add_argument("-o",
                        "--output-prefix",
                        type=str,
                        required=True,
                        help='output path (no extension, no ending .)')
    args = parser.parse_args()
    return args

def draw(in_file, out_file):
        f = open(in_file)
        # lines = f.readlines()

        colors = { 1 : 'b', 0 : 'g'}
        G = nx.Graph()
        for line in f:
            line = line.rstrip().split('\t')
            # for x in line[:-1]:
            #     print(float(x))
            # print(1)
            node_id = int(line[1])
            neighbor_count = int(line[2])
            read_count = int(line[3])

            left_avg = float(line[4])
            right_avg = float(line[5])
            hamming_avg = float(line[6])

            for neighbor in line[7].rstrip().split(','):
                if len(neighbor) < 1:
                    continue
                fields = neighbor.split('_')
                neighbor_id = int(fields[0])
                left = int(fields[1])
                right = int(fields[2])
                dist = int(fields[3])


                G.add_edge(node_id, neighbor_id, color=colors[dist],weight=(left+right)/5)

            plt.figure(figsize=(30,30))
            nx.draw(G)#, pos, edges=edges, edge_color=colors, width=weights)
            plt.savefig(out_file+"_basic.pdf")

            edges = G.edges()
            colors = [G[u][v]['color'] for u,v in edges]
            weights = [G[u][v]['weight'] for u,v in edges]

            plt.figure(figsize=(30,30))
            pos = nx.spring_layout(G)
            nx.draw(G, pos, edges=edges, edge_color=colors, width=weights, alpha=0.1, node_size=25)
            plt.savefig(out_file+"_spring.pdf")

            plt.figure(figsize=(18,18))
            pos = nx.spectral_layout(G)
            nx.draw(G, pos, edges=edges, edge_color=colors, width=weights)
            plt.savefig(out_file+"_spectral.pdf")

            # plt.figure(figsize=(18,18))
            # pos = nx.circular_layout(G)
            # nx.draw(G, pos, edges=edges, edge_color=colors, width=weights)
            # plt.savefig(out_file+"_circular.pdf")
            #
            # plt.figure(figsize=(18,18))
            # pos = nx.shell_layout(G)
            # nx.draw(G, pos, edges=edges, edge_color=colors, width=weights)
            # plt.savefig(out_file+"_shell.pdf")

            plt.figure(figsize=(30,30))
            pos = nx.random_layout(G)
            nx.draw(G, pos, edges=edges, edge_color=colors, width=weights)
            plt.savefig(out_file+"_random.pdf")

def main():
    args = parse_args()
    print("HI!")
    draw(args.input_debug_cluster, args.output_prefix)

if __name__ == "__main__":
    main()
