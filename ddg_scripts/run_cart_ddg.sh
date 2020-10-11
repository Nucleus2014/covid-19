#!/bin/bash
while (( $# > 1 ))
do
  case $1 in
    -tp) template_pdb="$2";;
    -ite) iterations="$2";;
    -part) partition="$2";;
    -mem) memory="$2";;
    *) break;
  esac; shift 2
done

if ! [ -z "${iterations}" ]
then
  iterations="-ddg:iterations "${iterations}
fi

if ! [ -z "${memory}" ]
then
  memory="--mem "${memory}
fi

for mut in `ls`
do
  cd ${mut:0:-4}
  slurmit.py --job ${mut:0:-4} --partition ${partition} --begin now ${memory} \
    --command "$ROSETTA3/bin/cartesian_ddg.linuxgccrelease -database $ROSETTA3_DB \
    -s ../${template_pdb} -ddg:mut_file ../${mut} ${iterations} -ddg::cartesian \
    -ddg::dump_pdbs true -ddg::bbnbrs 1 -fa_max_dis 9.0 -score:weights ref2015_cart_cst"
  sleep 0.05
  cd ..
done
exit