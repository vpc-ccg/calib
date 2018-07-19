#/bin/#!/usr/bin/env bash
 awk 'BEGIN{FS="\t"; OFS=" "} { if ($0 ~/^>/) {cid = substr($0, length(">Cluster 0"))} else {print cid,$2}}' $1 |\
  awk 'BEGIN{FS=" "; OFS="\t"} {n1 = substr($3,2,length($3)-4); print $1,"nid","rid",n1,"s1","q1","n2","s1","q1"}'
