#!/bin/bash
while (( $# > 1 ))
do
  case $1 in
    -tp) template_pdb="$2";;
    -ml) mutant_list="$2";;
    *) break;
  esac; shift 2
done

# Convert all ddg files in current directory into csv files
python3 ../../scripts/convert_ddg_to_csv.py

protein=${template_pdb%%"_"*}

if [[ "${mutant_list}" == *","* ]]
then
  cat ${protein}_1/${protein}_1.csv >> ${protein}.csv

  total_jobs=`ls ${protein}_*/${protein}_*.fingerprint | wc -l`
  for ((job_idx=2;job_idx<=total_jobs;job_idx++))
  do
    cat ${protein}_${job_idx}/${protein}_${job_idx}.csv >> ${protein}.csv
  done
fi

IFS=','
mutant_list=(${mutant_list[@]})
IFS=' '

for motif_idx in ${!mutant_list[@]}
do
  prefix=${mutant_list[$motif_idx]:0:-10}
  total_jobs=`ls ${prefix}_*/${prefix}_*.fingerprint | wc -l`
  for ((job_idx=1;job_idx<=total_jobs;job_idx++))
  do
    separated_csv=$(ls ${prefix}*_${job_idx}/${prefix}*_${job_idx}.csv)
    cat ${separated_csv} >> ${protein}.csv
  done
done

exit