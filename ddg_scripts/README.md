# commands

**Take monomer Orf6 as an example:**
1) Generate the pre-minimized .pdb file
> sbatch ../scripts/run_cart_min.sh Orf6_Reference.pdb (put cart2.script in the ../scripts directory)

2) Generate .mut and .fingerprint files
> ../scripts/generate_mutfile.sh -tp Orf6_cart.pdb -ml Orf6_GISAID.fasta.txt -ite 3 -wl 20 -part p_sdk94_1

-wl: Workload (decoys) for each slurm job.
-ite: Number of decoys for each variant, defalut 3.

3) Run cart ddg:
> ../scripts/run_cart_ddg.sh -tp Orf6_cart.pdb -ite 3 -part p_sdk94_1

4) Concatenate all *.ddg and *.fingerprint files:
> ../scripts/concatenate_saparated_ddg.sh -tp Orf6_cart.pdb -ml Orf6_GISAID.fasta.txt

**Take hetereodimer Nsp10-Nsp16 as an example:**
1) Generate the pre-minimized .pdb file
> sbatch ../scripts/run_cart_min.sh Nsp10-Nsp16_Reference.pdb (put cart2.script in the ../scripts directory)

2) Generate .mut and .fingerprint files
> ../scripts/generate_mutfile.sh -tp Nsp10-Nsp16_cart.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt -cut A,B -ite 3 -wl 20 -part p_sdk94_1

Nsp16_GISAID.fasta.txt corresponds to chain A while Nsp10_GISAID.fasta.txt corresponds to chain B.
Please keep the order of -ml and -cut unchanged in different calculations, 
since the order of the variants in 2 matched_0.fasta.txt files follows the order in the first fasta file.

3) Run cart ddg:
> ../scripts/run_cart_ddg.sh -tp Nsp10-Nsp16_cart.pdb -ite 3 -part p_sdk94_1

4) Concatenate all *.ddg and *.fingerprint files:
> ../scripts/concatenate_saparated_ddg.sh -tp Nsp10-Nsp16_cart.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt

**Take homodimer Nsp5 as an example:**
1) Generate the pre-minimized .pdb file
> sbatch ../scripts/run_cart_min.sh Nsp5_Reference.pdb (put cart2.script in the ../scripts directory)

2) Generate .mut and .fingerprint files
> ../scripts/generate_mutfile.sh -tp Nsp5_cart.pdb -ml Nsp5_GISAID.fasta.txt -cut A -dup A,B -ite 3 -wl 20 -part p_sdk94_1

-cut A: Nsp5_GISAID.fasta.txt is set to be aligned to chain A, so chain A now is the main chain.
-dup A,B: All duplicated chains corresponding to Nsp5_GISAID.fasta.txt.

3) Run cart ddg:
> ../scripts/run_cart_ddg.sh -tp Nsp5_cart.pdb -ite 3 -part p_sdk94_1

4) Concatenate all *.ddg and *.fingerprint files:
> ../scripts/concatenate_saparated_ddg.sh -tp Nsp5_cart.pdb -ml Nsp5_GISAID.fasta.txt
