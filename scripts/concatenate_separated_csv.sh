#!/bin/bash
while (( $# > 1 ))
do
 case $1 in
   -tp) template_pdb="$2";;
   -ml) mutant_list="$2";;
   *) break;
 esac; shift 2
done

prefix=${template_pdb%%"_"*}

cp ${prefix}_part1_mutants.csv ${prefix}_mutants.csv
cp ${prefix}_part1_substitutions.csv ${prefix}_substitutions.csv
total_jobs=`ls *_part*_mutants.csv | wc -l`
for i in  $(seq 2 ${total_jobs})
do
 tail -n +2 ${prefix}_part${i}_mutants.csv >> ${prefix}_mutants.csv
 tail -n +2 ${prefix}_part${i}_substitutions.csv >> ${prefix}_substitutions.csv
done

exit