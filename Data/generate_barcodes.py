import rstr

import random

random_seed = 42
number_of_barcodes = 1000
barcode_regex = r'[ACGT]{10}'
output_file_path = "barcodes"


barcodes = set()
while len(barcodes) < number_of_barcodes:
    barcodes.add(rstr.xeger(barcode_regex))

with open(output_file_path, 'w+') as output_file:
    for barcode in barcodes:
        print(barcode, file=output_file)
