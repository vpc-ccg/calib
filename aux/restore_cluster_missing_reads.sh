#/bin/#!/usr/bin/env bash
in_cluster=$1
fastq=$2
awk -v total_reads="$(wc -l $2 | awk '{print $1/4}')" '
    BEGIN{
        OFS="\t";
    }
    {
        read_ids[$3]++;
        if ($1 > max_cluster) {
            max_cluster = $1;
        }
        print;
    }
    END{
        for (i = 0; i < total_reads; i++) {
            if (!(i in read_ids)) {
                max_cluster++;
                print max_cluster,"nid",i,"n1","s1","q1","n2","s2","q2";
            }
        }
    }
    ' $1
