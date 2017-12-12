import random
import sys
import argparse
from itertools import islice


def parse_args():
    parser = argparse.ArgumentParser(
        description="Attaches barcodes from whitelist to ends of reads randomly")
    parser.add_argument("-b",
                        "--input-barcodes",
                        type=str,
                        required=True,
                        help="Barcodes txt file")
    parser.add_argument("-m",
                        "--input-molecules",
                        type=str,
                        required=True,
                        help="Molecules fasta file")
    parser.add_argument("-s",
                        "--random-seed",
                        type=int,
                        default=42,
                        help="NumPy random seed (default: 42)")
    parser.add_argument("-o",
                        "--output-barcoded-molecules",
                        type=str,
                        default=None,
                        help="Output file (default: stdout)")
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

            print(">{}_{}_{}".format(barcode_a, barcode_b, lines[0][1:]), file=output_file)
            print(barcode_a + lines[1] + barcode_b, file=output_file)


def main():
    args = parse_args()
    barcodes_file_path = args.input_barcodes
    molecules_file_path = args.input_molecules
    output_file_path = args.output_barcoded_molecules if args.output_barcoded_molecules else None
    seed = args.random_seed if args.random_seed else None
    attach_barcodes(barcodes_file_path, molecules_file_path, seed=seed, output_file_path=output_file_path)


if __name__ == "__main__":
    main()
