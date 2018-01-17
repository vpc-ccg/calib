awk 'BEGIN{c = 0} {if ($2 != c) {c++; printf"# "c"\n"} else {print "1\t2\t"$1}}' $1
