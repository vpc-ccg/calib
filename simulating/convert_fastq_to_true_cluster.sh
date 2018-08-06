awk 'BEGIN{OFS="\t"} NR%4 == 1 {printf $1"\t"; gsub(/:/,"\t"); print}' |
    awk 'BEGIN{OFS="\t"; FS="\t"; rid=0} {n1 = substr($1,2); cid=int(substr($2,2)); print cid,"nid",rid++,n1,"s1","q1","n2","s2","q2"}'
