import argparse
import rstr
import sys
import random


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generates a list of barcodes")
    parser.add_argument("-n",
                        "--num-of-barcodes",
                        type=int,
                        default=1000,
                        help="Number of barcodes (default: 1000)")
    parser.add_argument("-l",
                        "--len-of-one-end-barcode",
                        type=int,
                        default=10,
                        help="Length of one side of the barcode (default: 10)")
    parser.add_argument("-s",
                        "--random-seed",
                        type=int,
                        default=42,
                        help="NumPy random seed (default: 42)")
    parser.add_argument("-o",
                        "--output-barcodes",
                        type=str,
                        default=None,
                        help="Output barcodes file (default: stdout)")
    args = parser.parse_args()
    return args


def generate_barcodes(output_file_path=None,
                      random_seed=None,
                      number_of_barcodes=1000,
                      barcode_length=10):
    random.seed(random_seed)
    # barcode_regex = '[ACGT]{' +str(barcode_length) + '}'
    barcodes = set()
    while len(barcodes) < number_of_barcodes:
        barcode = "".join([random.sample(['A','C','G','T'], 1)[0] for _ in range(barcode_length)])
        barcodes.add(barcode)
    
    if output_file_path:
        output_file = open(output_file_path, 'w+')
    else:
        output_file = sys.stdout

    for barcode in barcodes:
        print(barcode, file=output_file)


def main():
    args = parse_args()
    generate_barcodes(output_file_path=args.output_barcodes, random_seed=args.random_seed,
                      number_of_barcodes=args.num_of_barcodes, barcode_length=args.len_of_one_end_barcode)


if __name__ == "__main__":
    main()
