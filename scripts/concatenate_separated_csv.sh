#!/bin/bash
while (( $# > 1 ))
do
  case $1 in
    -tp) template_pdb="$2";;
    -ml) mutant_list="$2";;
    *) break;
  esac; shift 2
done

protein=${template_pdb%%"_"*}

if [[ "${mutant_list}" == *","* ]]
then
  cat ${protein}_1_mutants.csv >> ${protein}_mutants.csv

  total_jobs=`ls ${protein}_*_mutants.csv | wc -l`
  for ((job_idx=2;job_idx<=total_jobs;job_idx++))
  do
    tail -n +2 ${protein}_${job_idx}_mutants.csv >> ${protein}_mutants.csv
  done
fi

IFS=','
mutant_list=(${mutant_list[@]})
IFS=' '

for motif_idx in ${!mutant_list[@]}
do
  prefix=${mutant_list[$motif_idx]:0:-10}

  if ! [[ -f ${protein}_mutants.csv ]]
  then
    echo $(head -n 1 ${prefix}_1_mutants.csv) >> ${protein}_mutants.csv
  fi
  
  total_jobs=`ls ${prefix}_*_mutants.csv | wc -l`
  for ((job_idx=1;job_idx<=total_jobs;job_idx++))
  do
    separated_csv=$(ls ${prefix}*_${job_idx}_mutants.csv)
    tail -n +2 ${separated_csv} >> ${protein}_mutants.csv
  done
done

exit
