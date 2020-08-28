# commands

**Take monomer Nsp1 as an example:**
1) Run the calculations. Generate site mutated PDB models and corresponding energies.

a) Three rounds of soft-rep repacking within 8.0 Å + fixed backbone soft-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp1_relaxed_INPUT.pdb -ml Nsp1_GISAID.fasta.txt -ite 3 -rep soft,soft -fix_bb true -rnd 3 -nbh 8.0 -wl 100 -part p_sdk94_1

b) Three rounds of soft-rep repacking within 8.0 Å + flexible backbone hard-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp1_relaxed_INPUT.pdb -ml Nsp1_GISAID.fasta.txt -ite 3 -rep soft,hard -rnd 3 -nbh 8.0 -wl 100 -part p_sdk94_1

c) Three rounds of hard-rep repacking within 8.0 Å + flexible backbone hard-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp1_relaxed_INPUT.pdb -ml Nsp1_GISAID.fasta.txt -ite 3 -rnd 3 -nbh 8.0 -wl 100 -part p_sdk94_1

d) Fast relax within 8.0 Å
> ../../scripts/run_on_amarel.sh -tp Nsp1_relaxed_INPUT.pdb -ml Nsp1_GISAID.fasta.txt -ite 3 -proto fastrelax -nbh 8.0 -wl 50 -part p_sdk94_1

2) Concatenate separated output .csv files into a single .csv file.

> ../../scripts/concatenate_separated_csv.sh -tp Nsp1_relaxed_INPUT.pdb -ml Nsp1_GISAID.fasta.txt

3) Analyze the output .csv file. See "analyze_results" folder.

**Take heterodimer Nsp10-Nsp16, which bears small organic molecule ligands as an example:**
1) Run the calculations. Generate site mutated PDB models and corresponding energies.

a) Three rounds of soft-rep repacking within 8.0 Å + fixed backbone soft-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp10-Nsp16_relaxed.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt -cut A,B -ite 3 -rep soft,soft -fix_bb true -rnd 3 -op true -nbh 8.0 -wl 100 -part p_sdk94_1

b) Three rounds of soft-rep repacking within 8.0 Å + flexible backbone hard-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp10-Nsp16_relaxed.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt -cut A,B -ite 3 -rep soft,hard -rnd 3 -op true -nbh 8.0 -wl 100 -part p_sdk94_1

c) Three rounds of hard-rep repacking within 8.0 Å + flexible backbone hard-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp10-Nsp16_relaxed.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt -cut A,B -ite 3 -rnd 3 -op true -nbh 8.0 -wl 100 -part p_sdk94_1

d) Fast relax within 8.0 Å
> ../../scripts/run_on_amarel.sh -tp Nsp10-Nsp16_relaxed.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt -cut A,B -ite 3 -proto fastrelax -op true -nbh 8.0 -wl 50 -part p_sdk94_1

2) Concatenate separated output .csv files into a single .csv file.

> ../../scripts/concatenate_separated_csv.sh -tp Nsp10-Nsp16_relaxed.pdb -ml Nsp16_GISAID.fasta.txt,Nsp10_GISAID.fasta.txt

3) Analyze the output .csv file. See "analyze_results" folder.

**Take heterotrimer Nsp12-Nsp7-Nsp8, which bears two stands of RNA as an example:**
1) Run the calculations. Generate site mutated PDB models and corresponding energies.

a) Three rounds of soft-rep repacking within 8.0 Å + fixed backbone soft-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp12-Nsp7-Nsp8_relaxed.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt -cut B,A,C -dup B,D -ite 3 -rep soft,soft -fix_bb true -rnd 3 -op true -nbh 8.0 -wl 100 -part p_sdk94_1

b) Three rounds of soft-rep repacking within 8.0 Å + flexible backbone hard-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp12-Nsp7-Nsp8_relaxed.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt -cut B,A,C -dup B,D -ite 3 -rep soft,hard -rnd 3 -op true -nbh 8.0 -wl 100 -part p_sdk94_1

c) Three rounds of hard-rep repacking within 8.0 Å + flexible backbone hard-rep minimization:
> ../../scripts/run_on_amarel.sh -tp Nsp12-Nsp7-Nsp8_relaxed.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt -cut B,A,C -dup B,D -ite 3 -rnd 3 -op true -nbh 8.0 -wl 100 -part p_sdk94_1

d) Fast relax within 8.0 Å
> ../../scripts/run_on_amarel.sh -tp Nsp12-Nsp7-Nsp8_relaxed.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt -cut B,A,C -dup B,D -ite 3 -proto fastrelax -op true -nbh 8.0 -wl 50 -part p_sdk94_1

2) Concatenate separated output .csv files into a single .csv file.

> ../../scripts/concatenate_separated_csv.sh -tp Nsp12-Nsp7-Nsp8_relaxed.pdb -ml Nsp8_GISAID.fasta.txt,Nsp12_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt

3) Analyze the output .csv file. See "analyze_results" folder.
