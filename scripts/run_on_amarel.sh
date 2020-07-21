#!/bin/bash
while (( $# > 1 ))
do
 case $1 in
   -tp) template_pdb="$2";;
   -ml) mutant_list="$2";;
   -sym) symmertry="$2";;
   -ex) ex_rotamers="$2";;
   -nbh) neighborhood_residue="$2";;
   -fr) fast_relax="$2";;
   -op) only_protein="$2";;
   -cut) cut_region_by_chains="$2";;
   -mem) membrane="$2";;
   -ind) ind_type="$2";;
   -debug) debugging_mode="$2";;
   -chk) check="$2";;
   -wl) workload="$2";;
   -part) partition="$2";;
   *) break;
 esac; shift 2
done

if ! [ -z "${symmertry}" ]
then
 symmertry="-sym ${symmertry}"
fi

if [ "${ex_rotamers}" == "true" ]
then
 ex_rotamers=""
else
 ex_rotamers="-no_ex"
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

if ! [ -z "${cut_region_by_chains}" ]
then
 cut_region_by_chains="-cut ${cut_region_by_chains}"
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

if [ "${check}" == "true" ]
then
 check="--check"
fi

if [ -z "${workload}" ]
then
 workload=15
fi

prefix=${template_pdb%%"_"*}

if [ -z "${cut_region_by_chains}" ]
then
 total_variants=$(expr `grep -o ">" ${mutant_list} | wc -l` - 1)
 if [ -z "${fast_relax}" ]
 then
  total_jobs=$((${total_variants} / ${workload} + 1))
 else
  total_jobs=$((${total_variants} * ${fast_relax} / ${workload} + 1))
  fast_relax="-fr ${fast_relax}"
 fi

 srun -J split_${mutant_list:0:-4} -p ${partition} -t 20:00 \
  python3 ../scripts/split_fasta.py -i ${mutant_list} -n ${total_jobs} \
  -t ${template_pdb} ${check}

 for job_idx in $(seq 1 ${total_jobs})
 do
  slurmit.py --job ${prefix}_job${i} --partition ${partition} --begin now \
   --command "python3 ../scripts/make_site_mutated_protein.py -t ${template_pdb} \
   -m ${mutant_list}_${template_pdb:0:-4}.${job_idx}.txt -rn ${prefix}_part${job_idx} \
   ${symmertry} ${ex_rotamers} ${neighborhood_residue} ${fast_relax} ${only_protein} \
   ${membrane} ${ind_type} ${debugging_mode}"
  sleep 0.1
 done
else
 if ! [ -z "${fast_relax}" ]
 then
  fast_relax="-fr ${fast_relax}"
 fi

 slurmit.py --job ${prefix} --partition ${partition} --begin now \
   --command "python3 ../scripts/make_site_mutated_protein.py -t ${template_pdb} \
   -m ${mutant_list} ${cut_region_by_chains} -rn ${prefix} ${symmertry} \
   ${ex_rotamers} ${neighborhood_residue} ${fast_relax} ${only_protein} \
   ${membrane} ${ind_type} ${debugging_mode}"
fi

exit
