#/bin/#!/usr/bin/env bash
awk 'BEGIN{FS="\t"} {print $3}' $1 |\
 awk 'BEGIN{FS=","; OFS="\t"} {for(i =1; i<=NF;i++) {$i-=1; print NR,"nid",$i,"n1","s1","q1","n2","s2","q2"}}'
