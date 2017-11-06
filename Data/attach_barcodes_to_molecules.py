import random
from itertools import islice

barcodes_file_path = "barcodes"
molecules_file_path = "molecules.fa"
output_file_path = "barcoded_molecules.fa"
random_seed = 42

with open(barcodes_file_path) as barcodes_file, open(molecules_file_path) as molecules_file, open(output_file_path, 'w+') as output_file:
    barcodes = [line.rstrip() for line in barcodes_file.readlines()]
    while(True):
        lines = [line.rstrip() for line in islice(molecules_file, 2)]
        if (len(lines)) != 2:
            break
        barcode_a = barcodes[random.randrange(0, len(barcodes))]
        barcode_b = barcodes[random.randrange(0, len(barcodes))]

        print(lines[0] + "_A:" + barcode_a + '_B:' + barcode_b, file=output_file)
        print(barcode_a + lines[1] + barcode_b, file=output_file)
