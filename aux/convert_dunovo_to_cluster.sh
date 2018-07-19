#/bin/#!/usr/bin/env bash
awk 'BEGIN{FS="\t"; OFS="\t"; cc = 0} {if (!($1 in cid)) {cid[$1]=cc++} print cid[$1],"nid","rid",$3,$4,$5,$6,$7$8}' $1 
