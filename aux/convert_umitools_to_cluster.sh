#/bin/#!/usr/bin/env bash
awk '
    BEGIN{
        OFS="\t"
    }
    NR > 1 {
        rid=int(substr($1, 1, match($1, "_") - 1));
        print $9,"nid",rid,"n1","s1","q1","n2","s2","q2";
    }
    ' $1
