import argparse
import rstr
import sys
import random

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generates whitelist of barcodes")
    parser.add_argument("-n", "--num-of-barcodes", type=int, help="Number of barcodes (default: 1000)",
                        default=1000)
    parser.add_argument("-o", "--output", type=str, help="Output barcodes file (default: stdout)")
    parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: no seed)", default=None)
    args = parser.parse_args()
    return args


def generate_barcodes(output_file_path=None,
                      random_seed=None,
                      number_of_barcodes=1000):
    random.seed(random_seed)
    barcode_regex = r'[ACGT]{10}'
    barcodes = set()
    while len(barcodes) < number_of_barcodes:
        barcodes.add(rstr.xeger(barcode_regex))

    if output_file_path:
        output_file = open(output_file_path, 'w+')
    else:
        output_file = sys.stdout

    for barcode in barcodes:
        print(barcode, file=output_file)


def main():
    args = parse_args()
    generate_barcodes(output_file_path=args.output, random_seed=args.random_seed,
                      number_of_barcodes=args.num_of_barcodes)


if __name__ == "__main__":
    main()
