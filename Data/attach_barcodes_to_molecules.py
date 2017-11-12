import random
import sys
import argparse
from itertools import islice


def parse_args():
    parser = argparse.ArgumentParser(
        description="Attaches barcodes from whitelist to ends of reads randomly")
    parser.add_argument("-i", "--input", type=str, help="Barcode whitelist path", required=True)
    parser.add_argument("-m", "--molecules", type=str, help="Molecules/amplicons fasta file", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output file (default: stdout)")
    parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: no seed)")
    args = parser.parse_args()
    return args


def attach_barcodes(barcodes_file_path, molecules_file_path, seed=None, output_file_path=None):
    with open(barcodes_file_path) as barcodes_file, open(molecules_file_path) as molecules_file:
        random.seed(seed)
        if output_file_path:
            output_file = open(output_file_path, 'w+')
        else:
            output_file = sys.stdout
        barcodes = [line.rstrip() for line in barcodes_file.readlines()]
        while True:
            lines = [line.rstrip() for line in islice(molecules_file, 2)]
            if (len(lines)) != 2:
                break
            barcode_a = barcodes[random.randrange(0, len(barcodes))]
            barcode_b = barcodes[random.randrange(0, len(barcodes))]

            print(">_A:" + barcode_a + '_B:' + barcode_b + '_' + lines[0][1:], file=output_file)
            print(barcode_a + lines[1] + barcode_b, file=output_file)


def main():
    args = parse_args()
    barcodes_file_path = args.input
    molecules_file_path = args.molecules
    output_file_path = args.output if args.output else None
    seed = args.random_seed if args.random_seed else None
    attach_barcodes(barcodes_file_path, molecules_file_path, seed=seed, output_file_path=output_file_path)


if __name__ == "__main__":
    main()
