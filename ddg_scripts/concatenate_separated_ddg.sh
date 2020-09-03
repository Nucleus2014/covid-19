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

  cp ${protein}_0.ddg2 ${protein}.ddg
  cp ${protein}_0.fingerprint ${protein}.fingerprint

  for motif_idx in ${!mutant_list[@]}
  do
    report_name_prefix=${protein}_${mutant_list[$motif_idx]:0:-10}_matched
    total_jobs=`ls ${report_name_prefix}_*.fingerprint | wc -l`
    for ((job_idx=1;job_idx<=total_jobs;job_idx++))
    do
      cat ${report_name_prefix}_${job_idx}.ddg2 >> ${protein}.ddg
      cat ${report_name_prefix}_${job_idx}.fingerprint >> ${protein}.fingerprint
    done
  done
else
  cp ${protein}_1.ddg2 ${protein}.ddg
  cp ${protein}_1.fingerprint ${protein}.fingerprint

  total_jobs=`ls ${protein}_*.fingerprint | wc -l`
  for ((job_idx=2;job_idx<=total_jobs;job_idx++))
  do
    cat ${protein}_${job_idx}.ddg2 >> ${protein}.ddg
    cat ${protein}_${job_idx}.fingerprint >> ${protein}.fingerprint
  done
fi
exit
