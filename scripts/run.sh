#!/bin/bash
while (( $# > 1 ))
do
 case $1 in
   -tp) template_pdb="$2";;
   -ml) mutant_list="$2";;
   -sym) symmertry="$2";;
   -ex) ex_rotamers="$2";;
   -nbh) neighborhood_residue="$2";;
   -op) only_protein="$2";;
   -cut) cut_region_by_chains="$2";;
   -mem) membrane="$2";;
   -ind) ind_type="$2";;
   -chk) check="$2";;
   -debug) debugging_mode="$2";;
   *) break;
 esac; shift 2
done

if ! [ -z "${symmertry}" ]
then
 symmertry="-sym "${symmertry}
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

if [ "${check}" == "true" ]
then
 check="--check"
fi

if [ "${debugging_mode}" == "true" ]
then
 debugging_mode="-debug"
else
 debugging_mode=""
fi

prefix=${template_pdb%%"_"*}

if [ -z "${cut_region_by_chains}" ] || [ ${#cut_region_by_chains} == 6 ]
then
 python ../scripts/split_fasta.py -i ${mutant_list} -t ${template_pdb} -n 12 ${check}

 for job_idx in {1..12}
 do
  python ../scripts/make_site_mutated_protein.py -t ${template_pdb} \
   -m ${mutant_list}_${template_pdb:0:-4}.${job_idx}.txt -rn ${prefix}_part${job_idx} \
   ${symmertry} ${ex_rotamers} ${neighborhood_residue} ${only_protein} \
   ${membrane} ${ind_type} ${debugging_mode} &
 done

 while ! [ ${completed_jobs} == 24 ]
 do
  completed_jobs=`ls *part*csv | wc -l`
  sleep 10
 done

 cp ${prefix}_part1_mutants.csv ${prefix}_mutants.csv
 cp ${prefix}_part1_substitutions.csv ${prefix}_substitutions.csv
 for i in {2..12}
 do
  tail -n +2 ${prefix}_part${i}_mutants.csv >> ${prefix}_mutants.csv
  tail -n +2 ${prefix}_part${i}_substitutions.csv >> ${prefix}_substitutions.csv
 done
else
 python ../scripts/make_site_mutated_protein.py -t ${template_pdb} \
  -m ${mutant_list} ${cut_region_by_chains} -rn ${prefix} ${symmertry} \
  ${ex_rotamers} ${neighborhood_residue} ${only_protein} ${membrane} \
  ${ind_type} ${debugging_mode}
fi

exit
