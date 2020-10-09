#!/bin/bash
#SBATCH --clusters=amarel
#SBATCH --partition=main
#SBATCH --job-name=cart
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm.%N.%j.log
#SBATCH --error=slurm.%N.%j.err
#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --begin=now
#SBATCH --open-mode=append

srun $ROSETTA3/bin/relax.linuxgccrelease \
  -database $ROSETTA3_DB \
  -s ${1} \
  -use_input_sc \
  -constrain_relax_to_start_coords \
  -ignore_unrecognized_res \
  -fa_max_dis 9.0 \
  -nstruct 1 \
  -score:weights ref2015_cart_cst \
  -relax:min_type lbfgs_armijo_nonmonotone

exit
