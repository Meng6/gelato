#!/bin/bash
###############################################################
# This script was generated by bashlog
# For more information, visit thomasrebele.org/projects/bashlog
###############################################################

export LC_ALL=C
mkdir -p data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384
rm -f data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/*
if type mawk > /dev/null; then awk="mawk"; else awk="awk"; fi
sort="sort "
check() { grep -- $1 <(sort --help) > /dev/null; }
check "--buffer-size" && sort="$sort --buffer-size=20% "
check "--parallel"    && sort="$sort --parallel=2 "

read_ntriples() { $awk -F" " '{ sub(" ", "\t"); sub(" ", "\t"); sub(/ \.$/, ""); print $0 }' "$@"; }
conv_ntriples() { $awk -F$'\t' '{ print $1 " " $2 " " $3 " ." }'; }



mkfifo data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/lock_mat0; ( $sort -t $'\t' -k 1 \
    <($awk -v FS=$'\t' ' 
          BEGIN { 
           out0c2_cond1["<advised_by>"] = "1"; 
          }
        
         (($2) in out0c2_cond1){ print $1 FS $3 } 
          ' ./data/raw/mgdb/bashlog/advised.tsv) > data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat0; mv data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/lock_mat0 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/done_mat0; cat data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/done_mat0 > /dev/null & exec 3> data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/done_mat0; exec 3>&-; ) & 

# plan
touch data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat1 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat2
$awk -v FS=$'\t' ' ($3 == "63244" && $2 == "<author>") { print $1 >> "data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat1" } 
 ($2 == "<author>") { print $1 FS $3 >> "data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat2" } 
  ' ./data/raw/mgdb/bashlog/dissertation.tsv 


mkfifo data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/lock_mat3; ( $sort -t $'\t' -k 1 \
    <(cat data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/lock_mat0 1>&2 2>/dev/null ;  \
        join -t $'\t' -1 1 -2 1 -o 1.2,2.2 \
        <($sort -t $'\t' -k 1 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat2) data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat0) > data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat3; mv data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/lock_mat3 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/done_mat3; cat data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/done_mat3 > /dev/null & exec 3> data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/done_mat3; exec 3>&-; ) & 

# plan
$sort -t $'\t' -k 1 -k 2 -u \
<(join -t $'\t' -1 1 -2 1 -o 1.1,2.2 \
    <($sort -t $'\t' -k 1 -u \
        <($awk -v FS=$'\t' '  { print $1} 
              ' \
            <($sort -t $'\t' -k 1 -k 2 -u \
                    <($awk -v FS=$'\t' '  { print $2 FS "63244"} 
                          ' \
                        <(cat data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/lock_mat0 1>&2 2>/dev/null ;  \
                            join -t $'\t' -1 1 -2 1 -o 1.1,1.2,2.1 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat0 \
                            <($sort -t $'\t' -k 1 -u data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat1))) \
                     | tee data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/full4 > data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/delta4
                while 
                
                $sort -t $'\t' -k 1 -k 2 -u \
                    <($awk -v FS=$'\t' '  { print $2 FS "63244"} 
                          ' \
                        <(cat data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/lock_mat3 1>&2 2>/dev/null ;  \
                            join -t $'\t' -1 1 -2 1 -o 1.1,1.2,2.1 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/mat3 \
                            <($sort -t $'\t' -k 1 -u \
                                <($awk -v FS=$'\t' '  { print $1} 
                                      ' data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/delta4)))) \
                     | comm -23 - data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/full4 > data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/new4;
                
                mv data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/new4 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/delta4 ; 
                $sort -u --merge -o data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/full4 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/full4 data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/delta4 ; 
                [ -s data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/delta4 ]; 
                do continue; done
                
                rm data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/delta4
                cat data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384/full4))) \
    <($sort -t $'\t' -k 1 \
        <($awk -v FS=$'\t' ' 
              BEGIN { 
               out0c2_cond1["<name>"] = "1"; 
              }
            
             (($2) in out0c2_cond1){ print $1 FS $3 } 
              ' ./data/raw/mgdb/bashlog/person.tsv)))

rm -rf data/query/mgdb/bashlog/unary_search_ancestors_for_63244_tmp8604320384
