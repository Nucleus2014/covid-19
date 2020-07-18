#!/bin/bash
while (( $# > 1 ))
do
 case $1 in
   -tp) template_pdb="$2";;
   -ml) mutant_list="$2";;
   -ex) ex_rotamers="$2";;
   -nbh) neighborhood_residue="$2";;
   -fr) fast_relax="$2";;
   -part) partition="$2";;
   *) break;
 esac; shift 2
done

if [ "${ex_rotamers}" == "true" ]
then
 ex_rotamers=""
else
 ex_rotamers="-no_ex"
fi
if ! [ -z "${neighborhood_residue}" ]
then
 neighborhood_residue="-nbh "${neighborhood_residue}
fi
if ! [ -z "${fast_relax}" ]
then
 fast_relax="-fr "${fast_relax}
fi
if ! [ -z "${partition}" ]
then
 partition="--partition "${partition}
fi

num=$(expr `grep -o ">" ${mutant_list} | wc -l` - 1)
python3 ../scripts/split_fasta.py -i ${mutant_list} -n ${num}

prefix=${template_pdb%%"_"*}

for i in $(seq 1 $num)
do
 mkdir ${prefix}"_variant"${i}
 mv ${mutant_list}"."${i}".txt" ${prefix}"_variant"${i}
 cd ${prefix}"_variant"${i}
 #slurmit.py --job ${prefix}_variant${i} ${partition} --begin now --command "python3 ../../scripts/relax_new_pdb.py -t ../"${template_pdb}" -m "${mutant_list}"."${i}".txt -rn "${prefix}"_variant"${i}" "${ex_rotamers}" "${neighborhood_residue}" "${fast_relax}
 python ../make_site_mutated_protease.py -t ../${template_pdb} -m ${mutant_list}.${i}.txt -rn ${prefix}_variant${i} ${ex_rotamers} ${neighborhood_residue} ${fast_relax}  
 cd ..
done
