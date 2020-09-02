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
  IFS=','
  mutant_list=(${mutant_list[@]})
  IFS=' '

  cp ${protein}_0_mutants.csv ${protein}_mutants.csv
  cp ${protein}_0_substitutions.csv ${protein}_substitutions.csv

  for motif_idx in ${!mutant_list[@]}
  do
    report_name_prefix=${protein}_${mutant_list[$motif_idx]:0:-10}
    total_jobs=`ls ${report_name_prefix}_*_mutants.csv | wc -l`
    for ((job_idx=1;job_idx<=total_jobs;job_idx++))
    do
      tail -n +2 ${report_name_prefix}_${job_idx}_mutants.csv >> ${protein}_mutants.csv
      tail -n +2 ${report_name_prefix}_${job_idx}_substitutions.csv >> ${protein}_substitutions.csv
    done
  done
else
  cp ${protein}_1_mutants.csv ${protein}_mutants.csv
  cp ${protein}_1_substitutions.csv ${protein}_substitutions.csv

  total_jobs=`ls ${protein}_*_mutants.csv | wc -l`
  for ((job_idx=2;job_idx<=total_jobs;job_idx++))
  do
    tail -n +2 ${protein}_${job_idx}_mutants.csv >> ${protein}_mutants.csv
    tail -n +2 ${protein}_${job_idx}_substitutions.csv >> ${protein}_substitutions.csv
  done
fi
exit
