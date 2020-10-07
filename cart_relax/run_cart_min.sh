#!/bin/bash
#SBATCH --clusters=amarel
#SBATCH --partition=p_sdk94_1
#SBATCH --job-name=orf7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3000
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm.%N.%j.log
#SBATCH --error=slurm.%N.%j.err
#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --begin=now
#SBATCH --open-mode=append

srun $ROSETTA3/bin/relax.linuxgccrelease \
  -s ${1} \
  -use_input_sc \
  -constrain_relax_to_start_coords \
  -ignore_unrecognized_res \
  -fa_max_dis 9.0 \
  -nstruct 20 \
  -relax:coord_constrain_sidechains  \
  -score:weights ref2015_cart \
  -relax:min_type lbfgs_armijo_nonmonotone \
  -relax:script ../scripts/cart2.script

exit
