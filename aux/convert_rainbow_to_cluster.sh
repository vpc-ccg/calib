awk 'BEGIN{c = "NA"} {if ($2 != c) {c=$2; printf"# "c"\n"} else {print "1\t2\t@"$1}}' $1
