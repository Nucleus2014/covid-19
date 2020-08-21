#!/bin/bash
#SBATCH --clusters=amarel
#SBATCH --partition=main
#SBATCH --job-name=ddgsym
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

srun /scratch/emd182/static_rosetta_builds/relax.static.linuxgccrelease \
  -s ${1} \
  -use_input_sc \
  -constrain_relax_to_start_coords \
  -ignore_unrecognized_res \
  -nstruct 20 \
  -relax:coord_constrain_sidechains  \
  -score:weights ref2015_cart \
  -relax:min_type lbfgs_armijo_nonmonotone \
  -relax:script /scratch/emd182/boot_camp/covid-19/scripts/cart2.script \
  -symmetry:symmetry_definition ${2} \
  -database /scratch/emd182/static_rosetta_builds/database \
  -out:suffix ${3}

exit
