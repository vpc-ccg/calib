import argparse

import matplotlib.pyplot as plt
import networkx as nx

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculates Rand Index for a clusters file")
    parser.add_argument("-i",
                        "--input-debug-cluster",
                        type=str,
                        required=True)
    parser.add_argument("-o",
                        "--output-prefix",
                        type=str,
                        required=False)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    f = open(args.input_debug_cluster)
    lines = f.readlines()

    G = nx.Graph()
