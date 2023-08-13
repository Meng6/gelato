#!/bin/bash
###############################################################
# This script was generated by bashlog
# For more information, visit thomasrebele.org/projects/bashlog
###############################################################

export LC_ALL=C
mkdir -p data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200
rm -f data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/*
if type mawk > /dev/null; then awk="mawk"; else awk="awk"; fi
sort="sort "
check() { grep -- $1 <(sort --help) > /dev/null; }
check "--buffer-size" && sort="$sort --buffer-size=20% "
check "--parallel"    && sort="$sort --parallel=2 "

read_ntriples() { $awk -F" " '{ sub(" ", "\t"); sub(" ", "\t"); sub(/ \.$/, ""); print $0 }' "$@"; }
conv_ntriples() { $awk -F$'\t' '{ print $1 " " $2 " " $3 " ." }'; }



mkfifo data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat0; ( $sort -t $'\t' -k 1 \
    <($awk -v FS=$'\t' ' 
          BEGIN { 
           out0c2_cond1["<advised_by>"] = "1"; 
          }
        
         (($2) in out0c2_cond1){ print $1 FS $3 } 
          ' ./data/raw/mgdb/bashlog/advised.tsv) > data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat0; mv data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat0 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat0; cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat0 > /dev/null & exec 3> data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat0; exec 3>&-; ) & 

mkfifo data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat1; ( $sort -t $'\t' -k 1 \
    <($awk -v FS=$'\t' ' 
          BEGIN { 
           out0c2_cond1["<name>"] = "1"; 
          }
        
         (($2) in out0c2_cond1){ print $1 FS $3 } 
          ' ./data/raw/mgdb/bashlog/person.tsv) > data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat1; mv data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat1 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat1; cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat1 > /dev/null & exec 3> data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat1; exec 3>&-; ) & 

# plan
touch data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat2 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat3
$awk -v FS=$'\t' ' ($3 == "63244" && $2 == "<author>") { print $1 >> "data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat2" } 
 ($2 == "<author>") { print $1 FS $3 >> "data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat3" } 
  ' ./data/raw/mgdb/bashlog/dissertation.tsv 


mkfifo data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat4; ( $sort -t $'\t' -k 1 \
    <(cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat0 1>&2 2>/dev/null ;  \
        join -t $'\t' -1 1 -2 1 -o 1.2,2.2 \
        <($sort -t $'\t' -k 1 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat3) data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat0) > data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat4; mv data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat4 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat4; cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat4 > /dev/null & exec 3> data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/done_mat4; exec 3>&-; ) & 

# plan
$sort -t $'\t' -k 1 -k 2 -k 3 -k 4 -u \
<(cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat1 1>&2 2>/dev/null ;  \
    join -t $'\t' -1 2 -2 1 -o 1.2,2.2,1.1,1.3 \
    <($sort -t $'\t' -k 2 \
        <(cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat1 1>&2 2>/dev/null ;  \
            join -t $'\t' -1 1 -2 1 -o 1.1,1.2,2.2 \
            <($sort -t $'\t' -k 1 \
                <($awk -v FS=$'\t' '  { print $1 FS $3} 
                      ' \
                    <($sort -t $'\t' -k 1 -k 2 -k 3 -u \
                            <($awk -v FS=$'\t' '  { print $2 FS "63244" FS "63244"} 
                                  ' \
                                <(cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat0 1>&2 2>/dev/null ;  \
                                    join -t $'\t' -1 1 -2 1 -o 1.1,1.2,2.1 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat0 \
                                    <($sort -t $'\t' -k 1 -u data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat2))) \
                             | tee data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/full5 > data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/delta5
                        while 
                        
                        $sort -t $'\t' -k 1 -k 2 -k 3 -u \
                            <($awk -v FS=$'\t' '  { print $2 FS "63244" FS $1} 
                                  ' \
                                <(cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/lock_mat4 1>&2 2>/dev/null ;  \
                                    join -t $'\t' -1 1 -2 1 -o 1.1,1.2,2.1 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat4 \
                                    <($sort -t $'\t' -k 1 -u \
                                        <($awk -v FS=$'\t' '  { print $1} 
                                              ' data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/delta5)))) \
                             | comm -23 - data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/full5 > data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/new5;
                        
                        mv data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/new5 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/delta5 ; 
                        $sort -u --merge -o data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/full5 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/full5 data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/delta5 ; 
                        [ -s data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/delta5 ]; 
                        do continue; done
                        
                        rm data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/delta5
                        cat data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/full5))) data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat1)) data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200/mat1)

rm -rf data/query/mgdb/bashlog/binary_search_ancestors_for_63244_tmp8664003200
