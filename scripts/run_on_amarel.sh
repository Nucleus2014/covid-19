#!/bin/bash
while (( $# > 1 ))
do
  case $1 in
    -tp) template_pdb="$2";;
    -ml) mutant_list="$2";;
    -cut) cut_region_by_chains="$2";;
    -ind) ind_type="$2";;
    -sym) symmertry="$2";;
    -mem) membrane="$2";;
    -rep) repulsive_type="$2";;
    -op) only_protein="$2";;
    -nbh) neighborhood_residue="$2";;
    -fix_bb) fix_backbone="$2";;
    -r) rounds="$2";;
    -fr) fast_relax="$2";;
    -debug) debugging_mode="$2";;
    -wl) workload="$2";;
    -part) partition="$2";;
    *) break;
  esac; shift 2
done

IFS=','
mutant_list=(${mutant_list[@]})
cut_region_by_chains=(${cut_region_by_chains[@]})
repulsive_type=(${repulsive_type[@]})
IFS=' '

if ! [ -z "${ind_type}" ]
then
  ind_type="-ind ${ind_type}"
fi

if ! [ -z "${symmertry}" ]
then
  symmertry="-sym ${symmertry}"
fi

if ! [ -z "${membrane}" ]
then
  membrane="-memb -mspan ${membrane}"
fi

if ! [ -z "${repulsive_type}" ]
then
  repulsive_type="-rep ${repulsive_type[@]}"
fi

if [ "${only_protein}" == "true" ]
then
  only_protein="-op"
else
  only_protein=""
fi

if ! [ -z "${neighborhood_residue}" ]
then
 neighborhood_residue="-nbh ${neighborhood_residue}"
fi

if [ "${fix_backbone}" == "true" ]
then
  fix_backbone="-fix_bb"
else
  fix_backbone=""
fi

if ! [ -z "${rounds}" ]
then
  rounds="-r "${rounds}
fi

if ! [ -z "${fast_relax}" ]
then
  fast_relax="-fr ${fast_relax}"
fi

if [ "${debugging_mode}" == "true" ]
then
  debugging_mode="-debug"
else
  debugging_mode=""
fi

if [ -z "${workload}" ]
then
  workload=15
fi

protein=${template_pdb%%"_"*}

if [ ${#mutant_list[@]} -gt 1 ]
then
  srun -J match_fasta -p ${partition} -t 20:00 \
    python3 ../scripts/match_fasta_replicates.py -i ${mutant_list[@]}

  for motif_idx in ${!mutant_list[@]}
  do
    mutant_list[$motif_idx]=${mutant_list[$motif_idx]:0:-10}"_matched_0.fasta.txt"
  done

  fastas=$( echo ${mutant_list[*]} )
  chains=$( echo ${cut_region_by_chains[*]} )
  slurmit.py --job ${protein}_0 --partition ${partition} --begin now \
    --command "python3 ../scripts/make_site_mutated_protein.py -t ${template_pdb} \
    -m ${fastas} -cut ${chains} -rn ${protein}_0 ${ind_type} \
    ${symmertry} ${membrane} ${repulsive_type} ${only_protein} ${neighborhood_residue} \
    ${fix_backbone} ${rounds} ${fast_relax} ${debugging_mode}"
  sleep 0.1

  for motif_idx in ${!mutant_list[@]}
  do
    mutant_list[$motif_idx]=${mutant_list[$motif_idx]:0:-20}
  done
else
  mutant_list[0]=${mutant_list[0]:0:-10}
fi

for motif_idx in ${!mutant_list[@]}
do
  total_variants=$(expr `grep -o ">" ${mutant_list[$motif_idx]}"_matched.fasta.txt" | wc -l` - 1)
  if [ -z "${fast_relax}" ]
  then
    total_jobs=$((${total_variants} / ${workload} + 1))
  else
    total_jobs=$((${total_variants} * ${fast_relax:4:} / ${workload} + 1))
  fi

  srun -J split_${mutant_list[$motif_idx]} -p ${partition} -t 20:00 \
    python3 ../scripts/split_fasta.py -i ${mutant_list[$motif_idx]}"_matched.fasta.txt" \
      -n ${total_jobs} -t ${template_pdb}
  rm ${mutant_list[$motif_idx]}"_matched.fasta.txt"

  if ! [ -z "${cut_region_by_chains}" ]
  then
    cut_region_by_chains[$motif_idx]="-cut ${cut_region_by_chains[$motif_idx]}"
  fi

  if [ ${#mutant_list[@]} -gt 1 ]
  then
    report_name_prefix=${protein}_${mutant_list[$motif_idx]}
  else
    report_name_prefix=${protein}
  fi

  for ((job_idx=1;job_idx<=total_jobs;job_idx++))
  do
    slurmit.py --job ${protein}_${job_idx} --partition ${partition} --begin now \
      --command "python3 ../scripts/make_site_mutated_protein.py -t ${template_pdb} \
      -m ${mutant_list[$motif_idx]}_matched_${job_idx}.fasta.txt \
      ${cut_region_by_chains[$motif_idx]} -rn ${report_name_prefix}_${job_idx} \
      ${ind_type} ${symmertry} ${membrane} ${repulsive_type} ${only_protein} \
      ${neighborhood_residue} ${fix_backbone} ${rounds} ${fast_relax} ${debugging_mode}"
    sleep 0.1
  done
done

exit
