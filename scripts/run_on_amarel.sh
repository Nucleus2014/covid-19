#!/bin/bash
while (( $# > 1 ))
do
  case $1 in
    -tp) template_pdb="$2";;
    -ml) mutant_list="$2";;
    -cut) cut_region_by_chains="$2";;
    -dup) duplicated_chains="$2";;
    -ind) ind_type="$2";;
    -sym) symmertry="$2";;
    -memb) membrane="$2";;
    -proto) protocol="$2";;
    -cart) cartesian="$2";;
    -rep) repulsive_type="$2";;
    -dis) fa_max_dis="$2";;
    -rep_wts) repulsive_weights="$2";;
    -nbh) neighborhood_residue="$2";;
    -fix_bb) fix_backbone="$2";;
    -op) only_protein="$2";;
    -rnd) rounds="$2";;
    -ite) iterations="$2";;
    -debug) debugging_mode="$2";;
    -wl) workload="$2";;
    -part) partition="$2";;
    -mem) memory="$2";;
    *) break;
  esac; shift 2
done

IFS=','
mutant_list=(${mutant_list[@]}) # convert mutant_list to array
cut_region_by_chains=(${cut_region_by_chains[@]})  # convert cut_region_by_chains to array
if ! [ -z "${repulsive_type}" ]
then
  repulsive_type="-rep "${repulsive_type[@]} # convert repulsive_type to string
fi
if ! [ -z "${repulsive_weights}" ]
then
  repulsive_weights="-rep_wts "${repulsive_weights[@]} # convert repulsive_weights to string
fi
if ! [ -z "${duplicated_chains}" ]
then
  duplicated_chains="-dup "${duplicated_chains[@]} # convert duplicated_chains to string
fi
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

if ! [ -z "${protocol}" ]
then
  protocol="-proto ${protocol}"
fi

if [ "${cartesian}" == "true" ]
then
  cartesian="-cart"
else
  cartesian=""
fi

if ! [ -z "${fa_max_dis}" ]
then
  fa_max_dis="-dis "${fa_max_dis}
fi

if ! [ -z "${iterations}" ]
then
  iterations="-ite ${iterations}"
fi

if ! [ -z "${rounds}" ]
then
  rounds="-rnd "${rounds}
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

if ! [ -z "${memory}" ]
then
  memory="--mem "${memory}
fi

protein=$(basename ${template_pdb%%"_"*})

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
  slurmit.py --job ${protein}_0 --partition ${partition} --begin now ${memory} \
    --command "python3 ../../scripts/make_site_mutated_protein.py -t ${template_pdb} \
    -m ${fastas} -cut ${chains} ${duplicated_chains} -rn ${protein}_0 ${ind_type} \
    ${symmertry} ${membrane} ${protocol} ${cartesian} ${fa_max_dis} ${repulsive_type} \
    ${repulsive_weights} ${neighborhood_residue} ${fix_backbone} ${only_protein} \
    ${rounds} ${iterations} ${debugging_mode} -no_cst_score"
  sleep 0.05

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
  if [ -z "${iterations}" ]
  then
    total_jobs=$(((${total_variants} + 1) / ${workload} + 1))
  else
    total_jobs=$(((${total_variants} + 1) * ${iterations:5} / ${workload} + 1))
  fi

  #srun -J split_${mutant_list[$motif_idx]} -p ${partition} -t 20:00 \
    python3 ../../scripts/split_fasta.py -i ${mutant_list[$motif_idx]}".fasta.txt" \
      -n ${total_jobs} -t ${template_pdb}

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
    slurmit.py --job ${report_name_prefix}_${job_idx} --partition ${partition} --begin now ${memory} \
      --command "python3 ../../scripts/make_site_mutated_protein.py -t ${template_pdb} \
      -m ${mutant_list[$motif_idx]}_${job_idx}.fasta.txt ${cut_region_by_chains[$motif_idx]} \
      ${duplicated_chains} -rn ${report_name_prefix}_${job_idx} ${ind_type} \
      ${symmertry} ${membrane} ${protocol} ${cartesian} ${fa_max_dis} ${repulsive_type} \
      ${repulsive_weights} ${neighborhood_residue} ${fix_backbone} ${only_protein} \
      ${rounds} ${iterations} ${debugging_mode} -no_cst_score"
    sleep 0.05
  done
done

exit
