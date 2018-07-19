#/bin/#!/usr/bin/env bash
awk 'BEGIN{FS="\t";OFS="\t"} {print $2,"nid",$1,"n1",$3,"q1","n2",$4,"q2"}' $1 
