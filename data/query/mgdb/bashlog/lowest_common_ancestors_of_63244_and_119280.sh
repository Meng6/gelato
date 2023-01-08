#!/bin/bash
###############################################################
# This script was generated by bashlog
# For more information, visit thomasrebele.org/projects/bashlog
###############################################################

export LC_ALL=C
mkdir -p tmp
rm -f tmp/*
if type mawk > /dev/null; then awk="mawk"; else awk="awk"; fi
sort="sort "
check() { grep -- $1 <(sort --help) > /dev/null; }
check "--buffer-size" && sort="$sort --buffer-size=34% "
check "--parallel"    && sort="$sort --parallel=2 "

read_ntriples() { $awk -F" " '{ sub(" ", "\t"); sub(" ", "\t"); sub(/ \.$/, ""); print $0 }' "$@"; }
conv_ntriples() { $awk -F$'\t' '{ print $1 " " $2 " " $3 " ." }'; }



mkfifo tmp/lock_mat0; ( $sort -t $'\t' -k 1 -u \
    <($awk -v FS=$'\t' ' 
          BEGIN { 
           out2_cond0c1["63244" FS "119280"] = "1"; 
          }
        
         (($1 FS $2) in out2_cond0c1){ print $3 } 
          ' ./data/query/mgdb/bashlog/interim_common_ancestors_of_63244_and_119280.tsv) > tmp/mat0; mv tmp/lock_mat0 tmp/done_mat0; cat tmp/done_mat0 > /dev/null & exec 3> tmp/done_mat0; exec 3>&-; ) & 

# plan
$sort -t $'\t' -k 1 -k 2 -u \
<(join  -v 1  -t $'\t' -1 1 -2 1 -o 1.1,1.2 \
    <($sort -t $'\t' -k 1 \
        <(cat tmp/lock_mat0 1>&2 2>/dev/null ;  \
            join -t $'\t' -1 1 -2 1 -o 1.1,2.2 tmp/mat0 \
            <($sort -t $'\t' -k 1 \
                <($awk -v FS=$'\t' ' 
                      BEGIN { 
                       out0c2_cond1["<name>"] = "1"; 
                      }
                    
                     (($2) in out0c2_cond1){ print $1 FS $3 } 
                      ' ./data/raw/mgdb/bashlog/person.tsv)))) \
    <($sort -t $'\t' -k 1 -u \
        <(cat tmp/lock_mat0 1>&2 2>/dev/null ;  \
            join -t $'\t' -1 1 -2 1 -o 1.2 \
            <($sort -t $'\t' -k 1 \
                <(cat tmp/lock_mat0 1>&2 2>/dev/null ;  \
                    join -t $'\t' -1 1 -2 1 -o 1.2,2.1 \
                    <($sort -t $'\t' -k 1 ./data/query/mgdb/bashlog/interim_advise.tsv) tmp/mat0)) tmp/mat0)))

 rm -f tmp/*