awk 'BEGIN{prev=""; n=1} {if ($1!=prev){prev =$1; print "# "n++} print "1\t2\t@"$3}' $1
