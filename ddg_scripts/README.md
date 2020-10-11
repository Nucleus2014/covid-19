# commands

**Take monomer Orf6 as an example:**
1) Generate the pre-minimized .pdb file and get the best decoy.

2) Generate .mut and .fingerprint files.
> ../../scripts/generate_mutfile.sh -tp Orf6_cart.pdb -ml Orf6_GISAID.fasta.txt -part p_sdk94_1

3) Run cart ddg.
> ../../scripts/run_cart_ddg.sh -tp Orf6_cart.pdb -part p_sdk94_1

4) Concatenate all output files.
> ../../scripts/concatenate_separated_ddg.sh -tp Orf6_cart.pdb -ml Orf6_GISAID.fasta.txt

5) Analyze the output .ddg file. See the analyze_results folder.

**Take hetereodimer Nsp10-Nsp16 as an example:**
1) Generate the pre-minimized .pdb file and get the best decoy.

2) Generate .mut and .fingerprint files.
> ../../scripts/generate_mutfile.sh -tp Nsp10-Nsp16_cart.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt -cut A,B -part p_sdk94_1

Nsp16_GISAID.fasta.txt corresponds to chain A while Nsp10_GISAID.fasta.txt corresponds to chain B.
Please keep the order of -ml and -cut unchanged in different calculations, 
since the order of the variants in 2 matched_0.fasta.txt files follows the order in the first fasta file.

3) Run cart ddg.
> ../../scripts/run_cart_ddg.sh -tp Nsp10-Nsp16_cart.pdb -part p_sdk94_1

4) Concatenate all output files.
> ../../scripts/concatenate_separated_ddg.sh -tp Nsp10-Nsp16_cart.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt

5) Analyze the output .ddg file. See the analyze_results folder.

**Take homodimer Nsp5 as an example:**
1) Generate the pre-minimized .pdb file and get the best decoy.

2) Generate .mut and .fingerprint files.
> ../../scripts/generate_mutfile.sh -tp Nsp5_cart.pdb -ml Nsp5_GISAID.fasta.txt -cut A -dup A,B -part p_sdk94_1

-cut A: Nsp5_GISAID.fasta.txt is set to be aligned to chain A, so chain A now is the main chain.
-dup A,B: All duplicated chains corresponding to Nsp5_GISAID.fasta.txt.

3) Run cart ddg.
> ../../scripts/run_cart_ddg.sh -tp Nsp5_cart.pdb -part p_sdk94_1

4) Concatenate all output files.
> ../../scripts/concatenate_separated_ddg.sh -tp Nsp5_cart.pdb -ml Nsp5_GISAID.fasta.txt

5) Analyze the output .ddg file. See the analyze_results folder.

**Take hetereotetramer Nsp12-Nsp7-Nsp8 as an example:**
1) Generate the pre-minimized .pdb file and get the best decoy.

2) Generate .mut and .fingerprint files.
> ../../scripts/generate_mutfile.sh -tp Nsp12-Nsp7-Nsp8_cart.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt -cut B,A,C -dup B,D -part p_sdk94_1

Nsp8_GISAID.fasta.txt corresponds to chain B and D, Nsp12_GISAID.fasta.txt corresponds to chain A while Nsp17_GISAID.fasta.txt corresponds to chain C.
Please keep the order of -ml and -cut unchanged in different calculations, 
since the order of the variants in 3 matched_0.fasta.txt files follows the order in the first fasta file.

3) Run cart ddg.
> ../../scripts/run_cart_ddg.sh -tp Nsp12-Nsp7-Nsp8_cart.pdb -part p_sdk94_1

4) Concatenate all output files.
> ../../scripts/concatenate_separated_ddg.sh -tp Nsp12-Nsp7-Nsp8_cart.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt

5) Analyze the output .ddg file. See the analyze_results folder.
