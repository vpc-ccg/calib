import argparse

import numpy
from pymsa.score import SumOfPairs
from pymsa.substitutionmatrix import SubstitutionMatrix

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prints sum of pair scores for MSA profiles")
    parser.add_argument("-f",
                        "--forward",
                        type=str,
                        required=True,
                        help="Input forward MSA profiles file")
    parser.add_argument("-r",
                        "--reverse",
                        type=str,
                        required=True,
                        help="Input reverse MSA profiles file")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Output file where different clusters sizes entropy stats are stored")
    args = parser.parse_args()
    return args

class PCR(SubstitutionMatrix):
    """ Costum substitution matrix for PCR duplication
    """
    def __init__(self, gap_penalty=-2):
        super(PCR, self).__init__(gap_penalty)
        self.distance_matrix = \
            {('A', 'A'): 0, ('A', 'C'): -1, ('A', 'G'): -1, ('A', 'T'): -1, ('A', 'N'): -1,
             ('C', 'A'): -1, ('C', 'C'): 0, ('C', 'G'): -1, ('C', 'T'): -1, ('C', 'N'): -1,
             ('G', 'A'): -1, ('G', 'C'): -1, ('G', 'G'): 0, ('G', 'T'): -1, ('G', 'N'): -1,
             ('T', 'A'): -1, ('T', 'C'): -1, ('T', 'G'): -1, ('T', 'T'): 0, ('T', 'N'): -1,
             ('N', 'A'): -1, ('N', 'C'): -1, ('N', 'G'): -1, ('N', 'T'): -1, ('N', 'N'): 0}
    def get_distance_matrix(self):
        return self.distance_matrix
    def get_distance(self, char1, char2):
        """ Returns the distance between two symbols

        :param char1:
        :param char2:
        :return: the distance value
        """
        if char1 is '-' and char2 is '-':
            result = 0
        elif char1 is '-' or char2 is '-':
            result = self.gap_penalty
        else:
            matrix = self.get_distance_matrix()
            if (char1, char2) in matrix:
                v = matrix[(char1, char2)]
            else:
                v = matrix[(char2, char1)]
            result = v
        return result


def main():
    args = parse_args()
    forward_file = open(args.forward)
    reverse_file = open(args.reverse)
    output_file = open(args.output, 'w+')

    scoring_obj = SumOfPairs(PCR())
    line_1 = forward_file.readline()
    line_2 = reverse_file.readline()
    sequences_1 = []
    sequences_2 = []
    while (len(line_1) > 0):
        if line_1[0] == '>':
            print(line_1.rstrip()[1:], file=output_file, end='\t')
            sequences_1 = []
            sequences_2 = []
        elif line_1[0] == '<':
            try:
                print(len(sequences_1), scoring_obj.compute(sequences_1)+scoring_obj.compute(sequences_2), sep='\t', file=output_file)
            except:
                for s in sequences_1:
                    print(s)
                print(len(sequences_1))
                for s in sequences_2:
                    print(s)
                print(len(sequences_2))
                exit()
        else:
            sequences_1.append(line_1.rstrip())
            sequences_2.append(line_2.rstrip())
        line_1 = forward_file.readline()
        line_2 = reverse_file.readline()

if __name__ == "__main__":
    main()
