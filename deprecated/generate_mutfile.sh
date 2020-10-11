#!/bin/bash
while (( $# > 1 ))
do
  case $1 in
    -tp) template_pdb="$2";;
    -ml) mutant_list="$2";;
    -cut) cut_region_by_chains="$2";;
    -dup) duplicated_chains="$2";;
    -ite) iterations="$2";;
    -wl) workload="$2";;
    -part) partition="$2";;
    *) break;
  esac; shift 2
done

IFS=','
mutant_list=(${mutant_list[@]}) # convert mutant_list to array
cut_region_by_chains=(${cut_region_by_chains[@]}) # convert cut_region_by_chains to array
if ! [ -z "${duplicated_chains}" ]
then
  duplicated_chains="-dup "${duplicated_chains[@]} # convert duplicated_chains to string
fi
IFS=' '

if [ -z "${iterations}" ]
then
  iterations=3
fi

if [ -z "${workload}" ]
then
  workload=15
fi

protein=${template_pdb%%"_"*}

if [ ${#mutant_list[@]} -gt 1 ]
then
  #srun -J match_fasta -p ${partition} -t 20:00 \
    python3 ../../scripts/match_fasta_replicates.py -i ${mutant_list[@]}

  for motif_idx in ${!mutant_list[@]}
  do
    mutant_list[$motif_idx]=${mutant_list[$motif_idx]:0:-10}"_matched_0.fasta.txt"
  done

  fastas=$( echo ${mutant_list[*]} )
  chains=$( echo ${cut_region_by_chains[*]} )
  slurmit.py --job ${protein}_0 --partition ${partition} --begin now \
    --command "python3 ../../scripts/generate_ddg_mutfile.py -t ${template_pdb} \
    -m ${fastas} -cut ${chains} ${duplicated_chains}"
  sleep 0.1

  for motif_idx in ${!mutant_list[@]}
  do
    mutant_list[$motif_idx]=${mutant_list[$motif_idx]:0:-20}"_matched"
  done
else
  mutant_list[0]=${mutant_list[0]:0:-10}
fi

for motif_idx in ${!mutant_list[@]}
do
  total_variants=$(expr `grep -o ">" ${mutant_list[$motif_idx]}".fasta.txt" | wc -l` - 1)
  total_jobs=$((${total_variants} * ${iterations} * 2 / ${workload} + 1))
  #srun -J split_${mutant_list[$motif_idx]} -p ${partition} -t 20:00 \
    python3 ../../scripts/split_fasta.py -i ${mutant_list[$motif_idx]}".fasta.txt" \
      -n ${total_jobs} -t ${template_pdb}

  if ! [ -z "${cut_region_by_chains}" ]
  then
    cut_region_by_chains[$motif_idx]="-cut ${cut_region_by_chains[$motif_idx]}"
  fi

  for ((job_idx=1;job_idx<=total_jobs;job_idx++))
  do
    slurmit.py --job ${protein}_${job_idx} --partition ${partition} --begin now \
      --command "python3 ../../scripts/generate_ddg_mutfile.py -t ${template_pdb} \
      -m ${mutant_list[$motif_idx]}_${job_idx}.fasta.txt \
      ${cut_region_by_chains[$motif_idx]} ${duplicated_chains}"
    sleep 0.1
  done
done

exit
