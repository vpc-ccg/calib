#!/bin/bash

samtools view -F 3852 -h "$1".sam |
    samtools sort -n -@ 16 -m 1G -T "$1".temp -o "$1".name_sorted.sam -O sam -;

awk '
{
    if (substr($0,1,1) == "@") {
        print > "'$1'.R1.sam";
        print > "'$1'.R2.sam";
    } else {
        # if ($2 < 256) {
            if (and($2, 64)) {
                print > "'$1'.name_sorted.R1.sam";
            } else {
                print > "'$1'.name_sorted.R2.sam";
            }
        # }
    }
}' "$1".name_sorted.sam;

paste \
 <(awk -f ./sam_to_pos.awk "$1".name_sorted.R1.sam) \
 <(awk -f ./sam_to_pos.awk "$1".name_sorted.R2.sam) |\
 awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}' > "$1".pos;

sort -S 2G -k 2,2 -k 5,5 -k 3,3n -k 6,6n "$1".pos > "$1".pos_sorted.pos;

awk '
function abs(v) {
    return v < 0 ? -v : v
}
BEGIN{
    f_readname = 1;
    f_chr1 = 2;
    f_chr2 = 5;
    f_pos1 = 3;
    f_pos2  = 6;
    f_seq1 = 4;
    f_seq2 = 7;
    prev_chr1 = -1;
    prev_chr2 = -1;
    prev_pos1 = -1;
    prev_pos2 = -1;
}
{
    if (prev_chr1 == $f_chr1 && prev_chr2 == $f_chr2 && abs($f_pos1 - prev_pos1) < 2 && abs($f_pos2 - prev_pos2) < 2) {
    } else {
        print "# "n++;
        # print "# "n-1, prev_chr1, $f_chr1, prev_chr2, $f_chr2, abs($f_pos1 - prev_pos1), abs($f_pos2 - prev_pos2)>  "awk.log";
        prev_chr1 = $f_chr1;
        prev_chr2 = $f_chr2;
        prev_pos1 = $f_pos1;
        prev_pos2 = $f_pos2;
    }
    # print prev_chr1, $f_chr1, prev_chr2, $f_chr2, abs(prev_pos1), abs(prev_pos2) > "awk.log";
    print $f_chr1":"$f_pos1":"$f_chr2":"$f_pos2"\t"r++"\t@"$f_readname"\t"$f_seq1"\tQ1\t@"$f_readname"\t"$f_seq2"\tQ2";
}' "$1".pos_sorted.pos > "$1".cluster
