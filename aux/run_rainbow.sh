rainbow=$1
forward_reads=$2
reverse_reads=$3
rainbow_mismatch=$4
rainbow_div=$5
rainbow_output=$6
if [ $rainbow_div = "false" ]; then
    $1 cluster \
        -1 $forward_reads \
        -2 $reverse_reads \
        -m $rainbow_mismatch \
        > $rainbow_output
else
    $1 cluster \
        -1 $forward_reads \
        -2 $reverse_reads \
        -m $rainbow_mismatch \
        | $1 div \
        > $rainbow_output
fi
