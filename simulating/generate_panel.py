import sys
import argparse
import random


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generates gene capture panel from GTF file and a gene list")
    parser.add_argument("-a",
                        "--gene-annotation",
                        type=str,
                        required=True,
                        help="Reference GTF file")
    parser.add_argument("-g",
                        "--gene-list",
                        type=str,
                        required=True,
                        help="Test file with gene name in each line")
    parser.add_argument("-n",
                        "--number-of-genes",
                        type=int,
                        default=30,
                        help="Number of random genes to select from the gene list to make the panel. (default: 35)")
    parser.add_argument("-s",
                        "--random-seed",
                        type=int,
                        default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("-o",
                        "--output-panel",
                        type=str,
                        default=None,
                        help="Output BED file of the generated random gene panel (default: stdout)")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    gene_names = list()
    for line in open(args.gene_list):
        gene_names.append(line.rstrip())

    random.seed(args.random_seed)
    random.shuffle(gene_names)
    gene_names = gene_names[0:args.number_of_genes]

    if args.output_panel:
        output_file = open(args.output_panel, 'w+')
    else:
        output_file = sys.stdout
    bed_records = list()
    gene_name_template = 'gene_name "{}";'
    for line in open(args.gene_annotation):
        if line[0] == "#":
            continue
        line = line.rstrip().split('\t')
        if line[2] != "exon":
            continue
        current_gene_name = "None"
        for gene_name in gene_names:
            if line[-1].find(gene_name_template.format(gene_name)) > 0:
                current_gene_name = gene_name
                break
        if current_gene_name == "None":
            continue
        current_chr = 'chr{}'.format(line[0].lstrip('chr'))
        current_start = str(int(line[3]))
        current_end = str(int(line[4]))
        bed_record = '\t'.join([current_chr, current_start, current_end, current_gene_name])
        bed_records.append(bed_record)
    for bed_record in bed_records:
        print(bed_record, file=output_file)

if __name__ == "__main__":
    main()
