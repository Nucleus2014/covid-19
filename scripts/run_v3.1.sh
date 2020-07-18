#!/bin/bash
while (( $# > 1 ))
do
 case $1 in
   -tp) template_pdb="$2";;
   -ml) mutant_list="$2";;
   -sym) symmertry="$2";;
   -ex) ex_rotamers="$2";;
   -nbh) neighborhood_residue="$2";;
   -ind) ind_type="$2";;
   -debug) debugging_mode="$2";;
   -chk) check="$2";;
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

python ../scripts/split_fasta_version3_align.py -i ${mutant_list} -n 12 -t ${template_pdb} ${check}

prefix=${template_pdb%%"_"*}

for job_idx in {1..12}
do
 python ../scripts/make_site_mutated_protease_v3.py -t ${template_pdb} -m ${mutant_list}_${template_pdb:0:-4}.${job_idx}.txt -rn ${prefix}_part${job_idx} ${symmertry} ${ex_rotamers} ${neighborhood_residue} ${ind_type} ${debugging_mode} &
 pids[${i}]=$!
done

for pid in ${pids[*]}
do
 wait $pid
done

for i in {1..12}
do
 rm ${mutant_list}.${i}.txt
done

mv ${prefix}_part1_mutants.csv ${prefix}_mutants.csv
mv ${prefix}_part1_substitutions.csv ${prefix}_substitutions.csv
for i in {2..12}
do
 tail -n +2 ${prefix}_part${i}_mutants.csv >> ${prefix}_mutants.csv
 rm ${prefix}_part${i}_mutants.csv
 tail -n +2 ${prefix}_part${i}_substitutions.csv >> ${prefix}_substitutions.csv
 rm ${prefix}_part${i}_substitutions.csv
done
