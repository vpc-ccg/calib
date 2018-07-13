# sort $1 -k1,1n -S 8G -t 8 | \
# awk -F "\t" 'BEGIN{cid = -1; rid = 0; OFS="\t"} {if ($1!=cid) {cid = $1; print "#",cid} print $2,rid++,$3,$4,$5,$6,$7,$8}'
awk -F "\t" '{if($1 in c) {c[$1] = c[$1]"\n"$0} else {c[$1]=$0}} END{for (k in c) {print "#\t"k; print c[k]}}' $1
