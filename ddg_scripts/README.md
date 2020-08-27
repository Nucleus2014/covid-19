# commands

**Take monomer Orf6 as an example:**
1) Generate the pre-minimized .pdb file and get the best decoy
> sbatch ../scripts/run_cart_min.sh Orf6_Reference.pdb (put cart2.script in the ../scripts directory)

> python ../scripts/get_the_lowest_cart_relax_decoy.py -sc score.sc

2) Generate .mut and .fingerprint files
> ../scripts/generate_mutfile.sh -tp Orf6_cart.pdb -ml Orf6_GISAID.fasta.txt -ite 3 -wl 30 -part p_sdk94_1

-wl: Workload (decoys) for each slurm job.
-ite: Number of decoys for each variant, defalut 3.

3) Run cart ddg:
> ../scripts/run_cart_ddg.sh -tp Orf6_cart.pdb -ite 3 -part p_sdk94_1

4) Concatenate all *.ddg and *.fingerprint files:
> ../scripts/concatenate_separated_ddg.sh -tp Orf6_cart.pdb -ml Orf6_GISAID.fasta.txt

5) Analyze the output .ddg file. See the analyze_results folder.

**Take hetereodimer Nsp10-Nsp16 as an example:**
1) Generate the pre-minimized .pdb file and get the best decoy
> sbatch ../scripts/run_cart_min.sh Nsp10-Nsp16_Reference.pdb (put cart2.script in the ../scripts directory)

> python ../scripts/get_the_lowest_cart_relax_decoy.py -sc score.sc

2) Generate .mut and .fingerprint files
> ../scripts/generate_mutfile.sh -tp Nsp10-Nsp16_cart.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt -cut A,B -ite 3 -wl 30 -part p_sdk94_1

Nsp16_GISAID.fasta.txt corresponds to chain A while Nsp10_GISAID.fasta.txt corresponds to chain B.
Please keep the order of -ml and -cut unchanged in different calculations, 
since the order of the variants in 2 matched_0.fasta.txt files follows the order in the first fasta file.

3) Run cart ddg:
> ../scripts/run_cart_ddg.sh -tp Nsp10-Nsp16_cart.pdb -ite 3 -part p_sdk94_1

4) Concatenate all *.ddg and *.fingerprint files:
> ../scripts/concatenate_separated_ddg.sh -tp Nsp10-Nsp16_cart.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt

5) Analyze the output .ddg file. See the analyze_results folder.

**Take homodimer Nsp5 as an example:**
1) Generate the pre-minimized .pdb file and get the best decoy
> sbatch ../scripts/run_cart_min_symm.sh Nsp5_Reference_INPUT.pdb (put cart2.script in the ../scripts directory)

> python ../scripts/get_the_lowest_cart_relax_decoy.py -sc score.sc

2) Generate .mut and .fingerprint files
> ../scripts/generate_mutfile.sh -tp Nsp5_cart.pdb -ml Nsp5_GISAID.fasta.txt -cut A -dup A,B -ite 3 -wl 30 -part p_sdk94_1

-cut A: Nsp5_GISAID.fasta.txt is set to be aligned to chain A, so chain A now is the main chain.
-dup A,B: All duplicated chains corresponding to Nsp5_GISAID.fasta.txt.

3) Run cart ddg:
> ../scripts/run_cart_ddg.sh -tp Nsp5_cart.pdb -ite 3 -part p_sdk94_1

4) Concatenate all *.ddg and *.fingerprint files:
> ../scripts/concatenate_separated_ddg.sh -tp Nsp5_cart.pdb -ml Nsp5_GISAID.fasta.txt

5) Analyze the output .ddg file. See the analyze_results folder.

**Take hetereotetramer Nsp12-Nsp7-Nsp8 as an example:**
1) Generate the pre-minimized .pdb file and get the best decoy
> sbatch ../scripts/run_cart_min.sh Nsp12-Nsp7-Nsp8_Reference.pdb (put cart2.script in the ../scripts directory)

> python ../scripts/get_the_lowest_cart_relax_decoy.py -sc score.sc

2) Generate .mut and .fingerprint files
> ../scripts/generate_mutfile.sh -tp Nsp12-Nsp7-Nsp8_cart.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt -cut B,A,C -dup B,D -ite 3 -wl 30 -part p_sdk94_1

Nsp8_GISAID.fasta.txt corresponds to chain B and D, Nsp12_GISAID.fasta.txt corresponds to chain A while Nsp17_GISAID.fasta.txt corresponds to chain C.
Please keep the order of -ml and -cut unchanged in different calculations, 
since the order of the variants in 3 matched_0.fasta.txt files follows the order in the first fasta file.

3) Run cart ddg:
> ../scripts/run_cart_ddg.sh -tp Nsp12-Nsp7-Nsp8_cart.pdb -ite 3 -part p_sdk94_1

4) Concatenate all *.ddg and *.fingerprint files:
> ../scripts/concatenate_separated_ddg.sh -tp Nsp12-Nsp7-Nsp8_cart.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt

5) Analyze the output .ddg file. See the analyze_results folder.
