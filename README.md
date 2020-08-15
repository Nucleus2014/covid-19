# Covid-19 project for klab
This project is operated by Sagar Khare[1][2], Stephen Burley[1], to discover mutations for each Covid-19 protein.

> [1] Institute for Quantitative Biomedicine, Rutgers University
>
> [2] Department of Chemistry, Rutgers University

Here is a update recording for the pipeline of generating mutations for covid-19 proteins, using GISAID database.

## Usage for scripts 

**Summarized by**: Changpeng Lu

### Help for Flags

#### 1. Must-include flags

* `-t`/`--template_pdb`: Input a starting PDB file for comparison and from which mutants will be generated.
* `-m`/`--mutants_list`: Input a fasta list file or multiple fasta files separated by space identifying the mutations.

#### 2. Choices of modeling protocols

* `-proto`/`--protocol`: Choose the modeling protocol, choices=['fastrelax', 'repack+min'], default='repack+min'.
* `-ite`/`--iterations`: Giving a flag of -ite, [ite] trajectories will be run for each variant and output the decoy with the lowest score.
* `-rep`/`--repulsive_type`: Choose score functions that will be used in repacking and minimization. Give two parameters in ['soft', 'hard'] separated by space.
* `-rnd`/`--rounds`: The rounds of repacking and minimization calculations being repeated in one trajectory.
* `-nbh`/`--neighborhood_residue`: Giving a flag of -nbh will also allow surrounding residues within [nbh] angstroms of the mutated residue to repack. Usually be the fastest way to do mutations.
* `-op`/`--only_protein`: Giving a flag of -op will prevent non-protein motifs from repacking, such as ligands and RNA.
* `-fix_bb`/`--backbone`: Giving a flag of -fix_bb will hold the backbone fixed in minimization.

#### 3. Specified flags required for different cases

* `-sym`/`--symmetry`: If the pose is symmetric, include a symdef file.
* `-memb`/`--membrane`: Declare if the protein is a membrane protein.
* `-mspan`/`--span_file`: If the pose is a membrane protein, include a spanfile.
* `-cut`/`--cut_region_by_chains`: if multiple fasta files input or only a part of regions to be mutated, regions are needed to be defined in the same order of fasta files order. example: "A C B"
* `-pmm`/`--is_pdb_index_match_mutant`: if pdb index matches mutant index, then set this flag to avoid pairwise alignment which costs much time and not so accurate for those alignments with many gaps.

#### 4. Other customized flags

* `-od`/`--out-dir`: Input a directory into which the homolog models will be saved.  If not specified, PDBs will be saved in the current directory.
* `-rn`/`--report_name`: Input a name for the substitutions report. If not specified, the report will be called substitutions_summary.csv and be saved in the current directory.
* `-cr`/`--catalytic_residues`: The catalytic residues of the enzyme. By default, no residues are so designated. If residues are specified, report will include whether substitutions interact with the catalytic residues.
* `-params`/`--params`: If a non-canonical residue/ligand is present, provide a params file.
* `-no_models`/`--make_models`: Giving a flag of -no_models will prevent PDB models from being generated. This will prevent the energetic calculations and not yield the substituted models, but will skip the most time-consuming steps.
* `-ind`/`--ind_type`: To show mutation residues indices in pdb or in pose order, choices are `pose` or `pdb`. Default: pdb.
* `-debug`/`--debugging_mode`: Giving a flag of -debug, it will print out the task operations on all residues, namely point mutations, repacking or keeping static.

### Examples for running the mutation script

*  Mutations on a single chain protein
   *  Example: nsp1
   *  `python make_site_mutated_protein.py -t nsp1_relaxed.pdb -m Nsp1_GISAID.fasta.txt -nbh 8.0 -rn Nsp1_fast`
*  Mutations on a Asymmetric protein with the same sequences for several chains
   * Example: 6VSB
   * `python make_site_mutated_protein.py -t 6VSB_relaxed.pdb -m S-protein_GISAID.fasta.txt S-protein_GISAID.fasta.txt S-protein_GISAID.fasta.txt -pmm -nbh 8.0 -rn 6VSB_fast -cut "A B C"`
*  Mutations on a asymmetric protein with different mutants
   *  Example: 7bv2
   *  `python make_site_mutated_protein.py -t 7bv2_relaxed.pdb -m Nsp12_GISAID.fasta.txt Nsp8_GISAID.fasta.txt Nsp7_GISAID.fasta.txt -pmm -nbh 8.0 -rn 7bv2_fast -cut "A B C"`
*  Mutations on only one chain of a asymmetric protein
   *  Example: 6m17
   *  `python make_site_mutated_protein.py -t 6m17_relaxed.pdb -m S-protein_GISAID.fasta.txt -nbh 8.0 -rn 6m17_fast -cut "E"`



## More information please go to the [log for developing pipeline and scripts](https://github.com/Nucleus2014/covid-19/blob/master/Log_covid-19_pipeline.md)


