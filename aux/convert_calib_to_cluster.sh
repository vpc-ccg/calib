sort $1 -k1,1n | \
awk -F "\t" 'BEGIN{cid = -1; rid = 0; OFS="\t"} {if ($1!=cid) {cid = $1; print "#",cid} print $2,rid++,$3,$4,$5,$6,$7,$8}' 
