#!/bin/bash
# 1 = rainbow
# 2 = forward_reads
# 3 = reverse_reads
# 4 = rainbow_mismatch
# 5 = rainbow_div
# 6 = rainbow_output
if [ $5 = "false" ]; then
    $1 cluster \
        -1 $2 \
        -2 $3 \
        -m $4 \
        > $6
else
    $1 cluster \
        -1 $2 \
        -2 $3 \
        -m $4 \
        | $1 div \
        > $6
fi
