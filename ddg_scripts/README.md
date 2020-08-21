# commands

**Take Nsp10-Nsp16 as a example:**
1) Generate the pre-minimized .pdb file
> sbatch ../scripts/run_cart_min.sh [pdb-name] (put cart2.script in the ../scripts directory)

2) Generate .mut and .fingerprint files
> ../scripts/generate_mutfile.sh -tp nsp10-nsp16_cart.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt -cut A,B -ite 3 -wl 20 -part p_sdk94_1

-wl: work load (decoys) for each slurm job -ite: number of decoys for each variant, defalut 3

3) Run cart ddg:
> ../scripts/run_cart_ddg.sh -tp nsp10-nsp16_cart.pdb -ite 3 -part p_sdk94_1

4) Concatenate all *.ddg and *.fingerprint files:
> ../scripts/concatenate_saparated_ddg.sh -tp nsp10-nsp16_cart.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt
