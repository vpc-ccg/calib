awk '{print $4}' $1 | \
awk -F ","  'BEGIN{c=0} {print "# "c; c++; for (i=1; i <= NF; i++) printf "1\t2\t@" $i-1"\n"}'
