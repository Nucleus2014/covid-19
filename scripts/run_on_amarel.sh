#!/bin/bash
while (( $# > 1 ))
do
  case $1 in
    -tp) template_pdb="$2";;
    -ml) mutant_list="$2";;
    -cut) cut_region_by_chains="$2";;
    -sym) symmertry="$2";;
    -nbh) neighborhood_residue="$2";;
    -fr) fast_relax="$2";;
    -op) only_protein="$2";;
    -mem) membrane="$2";;
    -ind) ind_type="$2";;
    -debug) debugging_mode="$2";;
    -wl) workload="$2";;
    -part) partition="$2";;
    *) break;
  esac; shift 2
done

IFS=','
read -ra mutant_list <<< ${mutant_list}
read -ra cut_region_by_chains <<< ${cut_region_by_chains}
IFS=' '

if ! [ -z "${symmertry}" ]
then
 symmertry="-sym ${symmertry}"
fi

if ! [ -z "${neighborhood_residue}" ]
then
 neighborhood_residue="-nbh ${neighborhood_residue}"
fi

if [ "${only_protein}" == "true" ]
then
 only_protein="-op"
else
 only_protein=""
fi

if ! [ -z "${membrane}" ]
then
 membrane="-memb -mspan ${membrane}"
fi

if ! [ -z "${ind_type}" ]
then
 ind_type="-ind ${ind_type}"
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

prefix=${template_pdb%%"_"*}

if [ ${#mutant_list[@]} > 1 ]
then
  srun -J match_fasta -p ${partition} -t 20:00 \
    python3 ../scripts/match_fasta_replicates.py -i ${mutant_list[@]}

  for motif_idx in ${!mutant_list[@]}
  do
    mutant_list[$motif_idx]=${mutant_list[$motif_idx]:0:-10}"_matched_0.fasta.txt"
  done

  slurmit.py --job ${prefix}_0 --partition ${partition} --begin now \
   --command "python3 ../scripts/make_site_mutated_protein.py -t ${template_pdb} \
   -m ${mutant_list[@]} -cut ${cut_region_by_chains[@]} -rn ${prefix}_0 \
   ${symmertry} ${neighborhood_residue} ${fast_relax} ${only_protein} \
   ${membrane} ${ind_type} ${debugging_mode}"
  sleep 0.1

  for motif_idx in ${!mutant_list[@]}
  do
    mutant_list[$motif_idx]=${mutant_list[$motif_idx]:0:-12}
  done
else
  mutant_list[0]=${mutant_list[0]:0:-10}
fi

for motif_idx in ${!mutant_list[@]}
do
  total_variants=$(expr `grep -o ">" ${mutant_list[$motif_idx]}".fasta.txt" | wc -l` - 1)
  if [ -z "${fast_relax}" ]
  then
    total_jobs=$((${total_variants} / ${workload} + 1))
  else
    total_jobs=$((${total_variants} * ${fast_relax} / ${workload} + 1))
    fast_relax="-fr ${fast_relax}"
  fi

  srun -J split_${mutant_list[$motif_idx]} -p ${partition} -t 20:00 \
    python3 ../scripts/split_fasta.py -i ${mutant_list[$motif_idx]}".fasta.txt" \
      -n ${total_jobs} -t ${template_pdb}

  if ! [ -z "${cut_region_by_chains}" ]
  then
    cut_region_by_chains[$motif_idx]="-cut ${cut_region_by_chains[$motif_idx]}"
  fi

  for job_idx in $(seq 1 ${total_jobs})
  do
    slurmit.py --job ${prefix}_${job_idx} --partition ${partition} --begin now \
      --command "python3 ../scripts/make_site_mutated_protein.py -t ${template_pdb} \
      -m ${mutant_list[$motif_idx]}"_"${job_idx}".fasta.txt" \
      ${cut_region_by_chains[$motif_idx]} -rn ${prefix}_${job_idx} \
      ${symmertry} ${neighborhood_residue} ${fast_relax} ${only_protein} \
      ${membrane} ${ind_type} ${debugging_mode}"
    sleep 0.1
  done
done

exit
