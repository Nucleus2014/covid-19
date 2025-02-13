# Covid-19 Pipeline Log

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
* `-m`/`--mutants_list`: Input a fasta list file or files identifying the mutations.

#### 2. Choices of mutation protocols

* `-nbh`/`--neighborhood_residue`: Giving a flag of -nbh will also allow surrounding residues within [nbh] angstroms of the mutated residue to repack. Usually be the fastest way to do mutations.
* `-fr`/`--fast_relax`: Giving a flag of -fr will employ fast relax protocol on whole protein instead of repacking and minimization, running for [fr] trajectories. It is a slowest way to do mutations.

#### 3. Specified flags required for different cases

* `-sym`/`--symmetry`: If the pose is symmetric, include a symdef file.
* `-memb`/`--membrane`: Declare if the protein is a membrane protein.
* `-mspan`/`--span_file`: If the pose is a membrane protein, include a spanfile.
* `-cut`/`--cut_region_by_chains`: if multiple fasta files input or only a part of regions to be mutated, regions are needed to be defined in the same order of fasta files order. example: "A,C,B"
* `-pmm`/`--is_pdb_index_match_mutant`: if pdb index matches mutant index, then set this flag to avoid pairwise alignment which costs much time and not so accurate for those alignments with many gaps.
* `-op`/`--only_protein`: Giving a flag of -op will prevent ligands and RNA motifs from repacking.

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
   * `python make_site_mutated_protein.py -t 6VSB_relaxed.pdb -m S-protein_GISAID.fasta.txt,S-protein_GISAID.fasta.txt,S-protein_GISAID.fasta.txt -pmm -nbh 8.0 -rn 6VSB_fast -cut "A,B,C"`
*  Mutations on a asymmetric protein with different mutants
   *  Example: 7bv2
   *  `python make_site_mutated_protein.py -t 7bv2_relaxed.pdb -m Nsp12_GISAID.fasta.txt,Nsp8_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt -pmm -nbh 8.0 -rn 7bv2_fast -cut "A,B,C"`
*  Mutations on only one chain of a asymmetric protein
   *  Example: 6m17
   *  `python make_site_mutated_protein.py -t 6m17_relaxed.pdb -m S-protein_GISAID.fasta.txt -nbh 8.0 -rn 6m17_fast -cut "E"`

## 2020-06-30 Main Script v1. release

**Released by**: 						  Joey Lubin

**Name of the main script**:  make_site_mutated_protease.py

First version is released by Joey Lubin, Rutgers University. The script mainly focus on detect where the substitutions between wild type sequence and the reference sequence in the FASTA, then make mutations on the PDB structure, repack and minimize to get the appropriate accepted low energy for the structure model.

**Functions of the code:**		Repacking + Minimization

*Note*: This repacking will be implemented to all residues for the protein.

**Command to run**:

```shell
python make_site_mutated_protease.py -t <respective pdb> -m <respective fasta file> -rn <respective nsp name>_fast -no_models
```

**More Instructions**:

* *In the case of Nsp5, I also added this flag to include checking against the catalytic residues: `-cr 41 145`*
* *For symmetry protein, please add `-sym <symfile.sym> flag*

* *Generating symmetry files uses the example command `rosetta/Rosetta/main/source/src/apps/public/symmetry/make_symmdef_file.pl -m CRYST -p 6lu7.pdb -r 4 > 6lu7.sym`*

* *Generating ligands params uses*

  *`python Rosetta/main/source/scripts/python/public/molfile_to_params.py -n T9J --chain B --root 1 -p T9J t9j.mol2`*

## 2020-06-30 Main Script Bugs fixed

**Fixed by**: Joey Lubin

Fix the bug that the mutations in the fasta file and csv file generated by your code do not match;

## 2020-06-30 Auxiliary patch v1. release 

**Released by**:								 Changpeng Lu

**Name of the patch**:					split_fasta.py

The script is to split the fasta file into N pieces that you could take for running separately on GANN server. N is suggested to be set the same as the number of nodes.

**Command to run**:

```shell
python split_fasta.py -i <fasta-file-name> -n <number of splits, default to 12>
```

**More instructions**:

* *Summarize all files into covid-19 folder for students to run, which makes things clearly*

* *With the help of Harry, installing pyrosetta could simply use `conda_env_create -f pyrosetta_conda/pyrosetta_environment.yml` on GANN servers. A more comprehensive tutorial that Elliott made is shown after*

* *Students are assigned to different hosts, which is shown below,*

  <img src="/Users/cplu/Desktop/Screen Shot 2020-07-10 at 5.15.46 AM.png" alt="Screen Shot 2020-07-10 at 5.15.46 AM" style="zoom:90%;" />

## 2020-07-01 Main Script Updates

**Updated by**: 						   Joey Lubin

**Name of the main script**:  make_site_mutated_protease.py

Add parallelization to the script, in case you don't want to split the FASTA files. If you want to use 12 processors, and this is the first processor, just add the following flag:

```shell
-part 12 1
```

This will split the list of fastas that the script iterates through into 12 segments and only iterate through the first one. The second integer should be a value between 1 and the first integer. Report files will have _1_of_12 appended to their names to prevent overwriting. I tested it, it will work fine with or without the parallelization flag.

## 2020-07-01 Tutorial v1. Release

**Released by**:					 Zhuofan shen

**Name of the tutorial**: 	tutorial.docx

Here is the step-by-step to run the scripts that Joey, Changpeng and I wrote. You may forward it to the students who have trouble running them. 

## 2020-07-03 Main Script v2. Release

**Released by**: 						 Zhuofan Shen

**Name of the main script**:  make_site_mutated_protease.py

I've updated several functions below in Joey's script so that it can run repacking+minimization or fast relax according to your choice.

* coord_constrain_pose(), according to Elliott's suggestion;

* make_point_mutant_task_factory(), allowing only residues nearby the mutations to repack;

* fast_relax_pose_with_muts(), fast relax a pose as an alternative to the "make_point_changes()" function, which do repacking+minimization;

* make_mutant_model(), for the reference pdb w/o any mutation, now this script will simply rescore the pose and do nothing else, because I can't set focus for NeighborhoodResidueSelector w/o any mutation.

**Command to run**:

```shell
python make_site_mutated_protease.py -t Nsp1_relaxed.pdb -m Nsp1_GISAID_fasta.txt -nbh 6.0 -fr 3
```

**More instructions**:

* *`-nbh` flag will set the distance to the mutations within which all residues are allowed to repack, leaving other residues static. If this flag is omitted, every residue would be able to repack. This flag takes effect in both fast relax and repacking+minimization operations, since both of them will use the task factory.*
* *`-fr` flag will tell the script to employ fast relax on the poses instead of repacking+minimization, and specify the numbers of decoys to run. If this flag is omitted, the script will simply do repacking+minimization on the poses.*
* *Plus, you can use the `-debug` flag to see which residues are mutated, repacked or stay static*
* *Besides, I also wrote a script to run fast relax jobs on Amarel. It splits a fasta file to each single variants, and then submit the relax job to slurm. To run this script, `cp -r my /home/zs251/covid-19` folder to your directory, `cd` to `nsp1` folder, and type: `bash run_relax_on_amarel.sh -tp Nsp1_relaxed.pdb -ml Nsp1_GISAID_fasta.txt -ex true -nbh 0.6 -fr 3 -part p_sdk94_1`*
* *The only problem is that "**PreventRepackingRLT**" operation doesn't work normally here*

## 2020-07-05 Main Script Bugs Fixed

**Fixed by**: Zhuofan Shen

It turns out that problem is from ***PyJobDistributor*** instead of ***PreventRepackingRLT***. To fix this issue, I let the job distributor output the relaxed decoys. So if you use the `--fast_relax` flag, every proteases variant originated from the same reference PDB should be placed into separate folders (or use abbreviated fasta file to run one variant).

## 2020-07-06 Tutorial v2. Release

**Released by**:					Zhuofan Shen

**Name of the tutorial**:   tutorial_v2.docx

## 2020-07-06 Auxiliary Updates

**Updated by**: Elliott Dolan

There seemed to be some confusion about getting pyrosetta on the gan servers. Joey and Ken spent some time putting together a yml file on them – which we thought wasn’t working.  I checked it out – it seems to work just fine.

Copy it from ken’s directory into your own home directory – 

- `/home/kdalenberg/pyrosetta_conda/pyrosetta_environment.yml`

Run the conda env create command below, pointing at the yml file.

- `conda env create -f pyrosetta_conda/pyrosetta_environment.yml` 

it takes about 10-20 minutes, but unpacks all the needed install stuff for the pyrosetta environment and it works perfectly as long as you activate the pyrosetta conda environment.

- `conda activate pyrosetta_env`

From here pyrosetta works. Open up a python terminal and `init()` and it works great.

## 2020-07-07 Auxiliary Script v2. Release

**Released by**: 				Changpeng Lu

**Name of the script **:   split_fasta_version2_for_bootcamp.py

In order to reach our goal more finely, I think both should be considered, clean fasta file by operating sequence alignments before running mutation, and detection of missing residues when doing mutation. This email is mainly focus on the former one. For now, using this newly version “split_fasta.py” script could solve mismatch between pdb and fasta file, by cleaning fasta files, assuming the pdb file is correct. 

So, I added a quality control part for my past “split_data.py” to check if the protein sequence in pdb file is different from the query one in fasta file. And I perform pairwise alignment among each pair of the pdb sequence with the sequence in fasta. The reason not to do multiple sequence alignment is that there is no need to install more alignment softwares such as ClustalW to generate MSA, and at the meantime, we could still receive adjusted version of fasta file directly from our past input fasta file that Joey made for us, only at the cost of 1-2 minutes computation time. 

Furthermore, this script combines the function of splitting fasta files and truncating sequences using alignments well. It could not only give you split files with truncated sequences without any more action, named after “.split_number.txt”,  but also the adjusted version of the fasta file with all mutant sequences, named after “.modified.txt”. Therefore, besides a little time consuming (around 1-2 min running 12 split files), all files we possibly need for the input can be generated.

**Command to run**:

```shell
python split_fasta_version2_for_bootcamp.py -i <original fasta file> -t <protein pdb file> -n <number of split, default 12> --check
```

**More Instructions**:

* *If input fasta file has different length from the one in pdb file, warning message will come out, but sequence alignments will not be performed unless set `--check` flag;*
* *If `-n` set to 1, then the output split file will only count to one file, and it will be the same content as the “modified” fasta file;*
* *Nsp13 as well as N-protein have been tested under this code and generated fasta files look nice. They are attached in this email so that you could use them directly without running the program again;*
* To make sure the problem for N-problem can be solved, I tried repacking + minimizing neighbor residues around 6 angstrom. It performs well and I attach the results “N-protein-results” to this email. 

* *One thing should be noticed by Zhuofan, when you operate one sequence per file, Joey’s script may have some technical output dimension mismatch problem. I tested on a file with a single reference sequence and his script got error. However, if I split all sequences into 12 files, it works very well except 11th fasta file. It is not because bugs in my code. The reason is that one of the mutant sequence which is `hCoV-19/Bosnia_and_Herzegovina/02-Banja_Luka`, lose some information on its tag. To make sure students could run it successfully, I offered the modified fasta files.*
* *Thanks to Zhuofan, he pointed out my script could not handle well when the total number of sequences could be divisible by the number of splits. It just because I added a if statement which is not necessary. Now the bug is fixed, and no additional actions should be addressed for users.*
* *Thanks to Elliott for suggesting me many commands and tips about doing alignments using pyrosetta modules, which makes my move faster.*

## 2020-07-09 Main Script v3. Release

**Released by**: 			  Changpeng Lu

**Name of the script**:  make_site_mutated_protease_v3.py

Here is the updated version3 of Joey’s script after I edited. Now, this script could handle original fasta files successfully, whatever the case is (listed below), as well as fasta files after alignment using my “split_fasta_version2_for_bootcamp.py”. The script will detect mismatches between pdb and fasta sequences automatically, correct them and show these mismatches in the output file.

**Command to run**:

```shell
python make_site_mutated_protease_v3.py -t <PDB filename> -m <FASTA file> -nbh 8.0 -rn <output name>_fast -ind pdb
```

**More Instructions**:

* *When considering symmetry protein, don’t forget to add `-sym <symmetry file.symm>`*
* *I used this script to generate N-protein output with pdb numbering, which has been attached. It seems good to me.*
* *For Nsp13 problem, I found that there is one deletion, so my past “split_fasta.py” fail to solve this problem. Therefore, Alan and Brea need to use the updated version of “make_site_mutated_protease_v3.py” to generate mutations. I generated outputs and they seem good to me.*

**Note**:

My updates on Joey’s script could handle such inputs when:

* Additions on N-terminal/C-terminal:

  For example, 

  > | POSE INDEX | **1**  | 2    | 3    | 4    | 5    | 6    |
  > | ---------- | ------ | ---- | ---- | ---- | ---- | ---- |
  > | PDB INDEX  | **24** | 25   | 26   | 27   | 28   | 29   |
  > | PDB seq    | **L**  | G    | T    | E    | P    | E    |
  > | FASTA seq  | **-**  | G    | T    | G    | P    | E    |

  Outputs:

  * Adjusted alignment

  > | PDB seq   | **L** | G    | T    | E    | P    | E    |
  > | --------- | ----- | ---- | ---- | ---- | ---- | ---- |
  > | FASTA seq | **L** | G    | T    | G    | P    | E    |

  * Mutation in CSV summary                 

    *  Two modes to choose: **pdb order** or **pose order (pose order will start from 1 and ignore gaps between regions)**

      **PDB mode:**

      ​		             **L24-; E27G**

      **POSE mode:**

  ​					                **L1-; E4G**

  * Output pdb filename for valid mutation

    ​		**PDB mode:**

    ​			                **N_hCoV-19_Australia_2043_2020_E27G.pdb**

    ​		**POSE mode;**

    ​			                **N_hCoV-19_Australia_2043_2020_E4G.pdb**

* Truncations on N-terminal/C-terminal:

  For example,

  > | POSE INDEX |       | 1    | 2    | 3     | 4    | 5    |
  > | ---------- | ----- | ---- | ---- | ----- | ---- | ---- |
  > | PDB INDEX  |       | 41   | 42   | 43    | 44   | 45   |
  > | PDB seq    | **-** | T    | A    | **Y** | P    | E    |
  > | FASTA seq  | **P** | T    | A    | **D** | P    | E    |

  Outputs,

  * Adjusted alignment

  > | POSE INDEX |      | 1    | 2    | 3     | 4    | 5    |
  > | ---------- | ---- | ---- | ---- | ----- | ---- | ---- |
  > | PDB INDEX  |      | 41   | 42   | 43    | 44   | 45   |
  > | PDB seq    |      | T    | A    | **Y** | P    | E    |
  > | FASTA seq  |      | T    | A    | **D** | P    | E    |
  * Mutation in CSV summary

    For truncation residues, I will only give a note that is two hyphens with index minus 1 between them.

    ​		**PDB mode:**

    ​			                **-40-; Y43D**

    ​		**POSE mode:**

    ​			                **-0-; Y3D**

  * Output filename for valid mutation

    ​		**PDB mode:**

    ​			                **N_hCoV-19_Australia_2043_2020_Y43D.pdb**

    ​		**POSE mode;**

    ​			                **N_hCoV-19_Australia_2043_2020_Y3D.pdb**

* Gaps

  For example,

  > | POSE INDEX | 1    | 2    | 3     | 4    | 5    | 6    |
  > | ---------- | ---- | ---- | ----- | ---- | ---- | ---- |
  > | PDB INDEX  | 24   | 25   | 26    | 27   | 28   | 29   |
  > | PDB seq    | L    | G    | **-** | -    | P    | E    |
  > | FASTA seq  | L    | G    | **T** | G    | P    | E    |

  Outputs,

  * Adjusted alignment: 

    keep it as before.

  * Mutation in CSV summary

    ​		**PDB mode:**

    ​			                **-26T; -27G**

    ​		**POSE mode:**

    ​			                **-3T; -4G**

  * Output filename for valid mutation

    No output pdbs

**The part of the code I edited is mainly on functions,**

* `compare_sequences(pdb_file, seq_2, ind_by)`, where pdb_file is the input of pdb file path, seq_2 is the sequence in FASTA, ind_by is record of index order mode (“pose” or “pdb”). It will return two lists of substitutions. One is for output with desired index order as well as those truncation/additions information, while the other one is for modelling.

* `analyze_mutant_protein(seqrecord, ref_pose, sf, ind_type, pdb_name, …)`. I only adjusted names of variables to keep logics of the script correct.

* I add one flag `-ind` | `—ind_type` to convert the mode of numbering. Choices are “pose” and “pdb”. The default is “pdb”.



Zhuofan made latest version of scripts, to run mutation, using

`../scripts/run_relax_on_amarel_v3_debug.sh -tp Nsp13_relaxed_INPUT.pdb -ml Nsp13_GISAID.fasta.txt -sym Nsp13_relaxed.symm -ex true -nbh 8.0 -fr 3 -ind pdb -debug true -chk true -wl 15 -part p_sdk94_1`

For students, the bash script is `run_v3.1.sh`

## 2020-07-09 Tutorial Updates

**Updated by**: Changpeng Lu

**Name of the tutorial**: Tutorial for N-protein mutation.html

This is the tutorial for running N-protein.

## 2020-07-10 Main Scripts Debugs Fixed

**Fixed by**: Zhuofan Shen

 I found a serious problem with the symmetric proteins when handling Luz's questions. After I modified the python script, it simply does a minimization for the wild type sequence, but the minimizer itself cannot effectively lower the energy of symmetric proteins. It seems that when Rosetta generates duplicated motifs for a symmetric protein, the starting structure is always bad. So I modified the script again, letting it do both repacking+minimization for the wild type sequence, and use the optimized wild type pose as the reference instead of the exact input wild type structure to calculate the energy change. Everything is perfect by now, but I would **not** suggest you use the **-nbh**" flag for big **symmetric** proteins. To begin with a bad initial structure, all residues in the wild type pose would be optimized, but the mutated variants would only be partly optimized. Comparing energy in this way makes no sense. I've updated all scripts in the Google Drive folder, hope you can redownload them to the server.

## 2020-07-10 Main Scripts Updated

**Updated By**: 				Changpeng Lu

**Name of the script**:   make_site_mutated_protease_v3.py

After talking with Sagar yesterday, we want to have CSV mutants output file to only show these that are substitutions (between the wild type and one of reference sequences), while CSV substitutions output file to show the list of these substitutions. Since the query sequence has already been aligned with reference sequences, the lengths of sequences from GISAID remain the same (I have checked all proteins that we have for now, nsp1, nsp2, nsp4,nsp5, nsp6, nsp9, nsp13, nsp14, nsp15, orf6, orf7a, N-protein.). Therefore, I aligned the PDB sequence with the reference sequence, one by one, and use index to find substitutions between the wild type and the reference.

**Outputs**:

* For CSV file to show mutants, the example table is shown below.

| Date_first | Date_last | count | Id_1                          | Id_2                          | Location_1 | Location_2 | tag  | substitutions         | Is_in_PDB         | multiple | conservative     | Energy_change | rmsd       | layer | coevolution |
| ---------- | --------- | ----- | ----------------------------- | ----------------------------- | ---------- | ---------- | ---- | --------------------- | ----------------- | -------- | ---------------- | ------------- | ---------- | ----- | ----------- |
| 3/25/20    | 3/25/20   | 1     | hCoV-19/Australia/VIC597/2020 | hCoV-19/Australia/VIC597/2020 | Australia  | Australia  | N    | D0Y;R203K;G204R;C999A | False;False;False | FALSE    | False;True;False | -0.2985954    | 0.04345622 |       | FALSE       |

1. “**is_in_PDB**” is a new column I added. It will show whether there is a mutation we could make in PDB. If False, it means PDB protein has a gap for this position, or N/C terminal region for the reference is larger than that in PDB;

2. Substitutions numbering logic:

   If the number in the middle is 0, such as “D0Y”, it means this substitution happens beyond the N-terminal of PDB protein;

   If the number in the middle is 999, such as “C999A”, it means this substitution happens beyond the C-terminal of PDB protein;

* For CSV file to show substitutions, the example table is shown below.

| site | native | mutant | Is_in_PDB |
| ---- | ------ | ------ | --------- |
| 0    | G      | A      | FALSE     |
| 62   | E      | G      | TRUE      |
| 999  | T      | N      |           |

The metric is similar to mutants CSV file.

**Note**:

* For now, if substitutions are beyond the N/C-terminal of PDB protein, the site will be labeled as 0 or 999 for all substitutions. Do you want to know the exact numbering records for those substitutions? Maybe count those beyond the N-terminal with “-1”, “-2”, “-3”? And for those beyond the C-terminal, just numbering after the final residue index in PDB?
* I saw Zhuofan has updated the script to handle symmetry proteins better. I am not sure if I need to add it as well, so I didn’t consider his updates for now. I also didn’t generate results for symmetry proteins. Do you want me to include Zhuofan’s update as well and generate mutations for those symmetric proteins?

## 2020-07-10 Main Script v4 Release

**Released by**: 				Changpeng Lu

**Name of the script**: 	make_site_mutated_protease_v4.py

As sagar mentioned early in the morning, which is shown below,

>  "Could we have the number inferred from alignment with PDB instead of 0 or 999 (i.e. if PDB starts at 89, residue in Fasta will be that is currently “0” will be 88, 87 etc.). Similarly for C-term."

Now the script will number those substitutions beyond both N and C terminals. Now it could fix renumbering problems in a fine way. To distinguish it with the former changes, I rename it as version 4 and please use this version for now to run mutations.

## 2020-07-10 Main Script v4 debug fixed

**Fixed by**: Changpeng Lu

Here is the updated version of mutation version 4. It could solve the problem with recording energies for non-structured residues too.

## 2020-07-11 Auxiliary Script v3 debug fixed

**Reported by**: Lingjun Xie

With Zhuofan’s help, now I can run the scripts on Orf8, PLpro and nsp3e, but for PLpro, I can’t get the part 9, and the final csv file only contains the 12th part. For nsp3e, I can get all the pieces, but again, the final csv file only contains the 12th part.

For the Orf7b, I keep getting this error:

```bash
Traceback (most recent call last):

 File "../scripts/split_fasta_version3_align.py", line 125, in <module>

  gp.write(que_out[(batch + 1) * 60 :] + "\n")     

NameError: name 'batch' is not defined
```

The command I use is: `../scripts/run_v4.sh -tp Orf7b_relaxed.pdb -ml Orf7b-GISAID.fasta.txt -ex true -ind pdb`



**Fixed by**: Changpeng Lu

Split script was not able to handle small protein. Now I fixed this problem and try it on Orf7b. The splited files seem good. 

## 2020-07-12 Main Script v4 debug fixed

**Reported by**: Vidur Sarma

The calculations for nsp3_nsp1a ran alright and the mutant files were generated. But I’ve run into the following error for nsp3_nsp3b and nsp3_UNK. The output stops at the first sequence, so I’m not sure whether this is a sequence alignment issue or something else. It would great if you could look into it and help out.

```
Traceback (most recent call last):

 File "../scripts/make_site_mutated_protease_v4.py", line 1031, in <module>

  main(args)

 File "../scripts/make_site_mutated_protease_v4.py", line 975, in main

  ind_type = args.ind_type, pdb_name = args.template_pdb)

 File "../scripts/make_site_mutated_protease_v4.py", line 858, in analyze_mutant_protein

  substitutions, new_subs = compare_sequences(pdb_name, seqrecord.seq, query, ind_type)

 File "../scripts/make_site_mutated_protease_v4.py", line 204, in compare_sequences

  start_ind = int(tmp[5]) # record the start index of residues

ValueError: invalid literal for int() with base 10: '24.972'
```

**Fixed by**: Changpeng Lu

I record the starting index for pdb sequence by splitting rows in the pdb files using space, however for nsp3b and nsp_UNK, the residue indices are followed right after chain name. Also, I found that “\n” is enough for separating two aligned sequences in the alignment results and the former one I used (“\n “) will get error when the tag of reference sequence doesn’t start from a space “ “. I have updated the version attached and I also uploaded it to two drives folder. I checked each test sample from nsp3_nsp3b and nsp3_UNK, and they work well.

## 2020-07-13 Main Script v4 debug fixed

**Reported by**: Lingjun and Zhuofan

PL_Pro isn't working.

```bash
Traceback (most recent call last):
  File "make_site_mutated_protease_v4.py", line 1032, in <module>
    main(args)
  File "make_site_mutated_protease_v4.py", line 976, in main
    ind_type = args.ind_type, pdb_name = args.template_pdb)
  File "make_site_mutated_protease_v4.py", line 854, in analyze_mutant_protein
    mut_tags = read_name_tag(seqrecord.id)
  File "make_site_mutated_protease_v4.py", line 182, in read_name_tag
    tag_dict['tag'] = breakup_id[5]
IndexError: list index out of range
```

**Fixed by**: Changpeng Lu

The reason why PLpro isn’t working is because the original GISAID file has some flaws in it. hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536 doesn’t have tag in its label. Need to be corrected by hand.

If it happens again, please  first look at the tag of the sequence which is shown before this error, which is

` ['2020-02-01', '2020-02-01', 'Count=1', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/2020', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/20']`

If the tag is right, this list above should have 5 elements, which should be `['2020-02-01', '2020-02-01', 'Count=1', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/2020', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/20', ‘NSP3_PL_PRO’]`

However in this case, it could not be processed correctly. So you need to locate this sequence in the FASTA file, which is

`2020-02-01|2020-02-01|Count=1|hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/2020|hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/20`

So if you add **“|NSP3_PL_PRO”** after the tag of sequence, the problem should be fixed.

## 2020-07-13 Main Script v4 debug fixed

**Reported by**: Sagar Khare

Just find out the problem is that mutation beyond pdb cannot be shown in the results. As long as you have it in mutants.csv it is fine to not have it in substitutiins.csv



**Fixed by **: Changpeng Lu

There is a small mistake I made so that the script will not output those mutations out of pdb sequence. Sorry for that! I uploaded the updated version already and tested on nsp3e attached. They look good to me.



## 2020-07-14 Auxiliary Script v3 debug fixed

**Fixed by**: Changpeng Lu

The problem is because splited fasta files generated from split script truncate reference sequences to match pdb sequences, so when mutation script analyzes alignment between pdb and reference sequences, those regions beyond the pdb has already been truncated and can never be detected. So I modified a little bit about split script and you need to use another command below:

`python split_fasta_version3_align.py -i NSP3_NSP3B-gisaid-20200625.fasta.txt -t Nsp3_Nsp3b_1025_1230_relaxed_0_0_INPUT.pdb`

I checked on one of nsp3b splited fasta files and now it works. Modified scripts are attached. I also uploaded on the drive.

## 2020-07-15 Main and Auxiliary Script New Version Released

**Released by**: Changpeng Lu

**Name of the script**: 	make_site_mutated_protease_v5.py; 

​											split_fasta_version4.py

**A general view is, now two scripts are updated to the new version. For split script, is “split_fasta_version4.py”; for mutation script, is “make_site_mutated_protease_v5.py”. Please wait for Zhuofan’s new version bash script.**

1. If you found there are     no outputs in **substitution** column when "**FALSE**" **is_in_PDB** column **at all**, that means the updated mutation script has not been used in your process. 

For now, it is mainly because Zhuofan didn't update his bash script to handle the updated mutation script. Please wait for a little bit for his updates;

For those who want to give it a shot without executed bash script, you could use command **without “—check”** that I mentioned in the group meeting to run the “split_fasta_version4.py” to split the fasta and run the mutation script as usual;

2. If you found some have outputs in the **substitution** column when "**FALSE**" in **is_in_PDB** column, but it cannot match your alignment analysis using BLAST ([https://blast.ncbi.nlm.nih.gov/Blast.cgi](https://slack-redir.net/link?url=https%3A%2F%2Fblast.ncbi.nlm.nih.gov%2FBlast.cgi)), that is mainly my fault on split_fasta script. The main reason is that the split_fasta script assigned a wrong series of reference sequences to the part 12 fasta file. Now it has been solved, and split fasta file is updated to the **version4**, but Zhuofan needs some time to correct his bash     script. For now, you could check your other proteins' outputs to see if they have the same issue. The way to check is to see if substitutions in     mutants.csv for the last several reference proteins are the same as those for the first several reference proteins in the same file. If it happens,  you could re-split fasta file using the command I just mentioned in point1 and see if things could be fixed.

 

3. Some reported the error below:

```bash
['2020-02-01', '2020-02-01', 'Count=1', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/2020', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/20']

Traceback (most recent call last):

 File "make_site_mutated_protease_v4.py", line 1032, in <module>

  main(args)

 File "make_site_mutated_protease_v4.py", line 976, in main

  ind_type = args.ind_type, pdb_name = args.template_pdb)

 File "make_site_mutated_protease_v4.py", line 854, in analyze_mutant_protein

  mut_tags = read_name_tag(seqrecord.id)

 File "make_site_mutated_protease_v4.py", line 182, in read_name_tag

  tag_dict['tag'] = breakup_id[5]

IndexError: list index out of range
```



I have mentioned this error in the former email. That is basically because hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536 doesn’t have tag in its label in the original FASTA file. I don’t think our script needs to handle this problem because it is very rare and could be quickly fixed in hand.


 The way to solve this, is to first look at the tag of the sequence which is shown before this error, ['2020-02-01', '2020-02-01', 'Count=1', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/2020', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/20']If the tag is right, this list above should have 5 elements, which should be ['2020-02-01', '2020-02-01', 'Count=1', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/2020', 'hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/20', **‘NSP3_PL_PRO’**]

 

When looking at this sequence in the FASTA file,
 *>2020-02-01|2020-02-01|Count=1|hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/2020|hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/20*


 It could not be processed correctly. **The good example is shown below,****
** >2020-02-01|2020-02-01|Count=1|hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/2020|hCoV-19/Germany/BavPat15-ChVir1484-ChVir1536/20**|NSP3_PL_PRO**

 

So if you add **“|NSP3_PL_PRO”** after the tag of sequence, the problem should be fixed. If not, please send us the message. I will take a look.

 

**TIP:** If you have some bugs, maybe it is already been posted by others. I tracked all these bugs in **Log_covid-19_pipeline** so you could quickly search if your bugs has already been fixed by finding **“debug”** word in the file. That is why I created this file. You could look at it in the drive: https://drive.google.com/drive/u/1/folders/1Rw5t6VrcWvrum-fF_CUABdBx7lU-8Svw. And report your bugs on it for others to follow easily. 



## 2020-07-15 Execution Script v4.1 Released

**Released by**: Zhuofan

Hi everyone, thanks to Changpeng's effort, the split_fasta_v4.py script works well now. She also remided me that I shouldn't use the --check flag while running the spliting script. So I updated both the shell script to run4.1.sh. You can download both the scripts from google drive.

# 2020-07-16 Main Script Debugs Summary

**Reported by**: Vidur Sarma

As we had discussed in the group meeting, the spike-Ace2 complex has a whole chain more than the fasta sequence i.e. the input structure has about 375 residues more than in the fasta sequence. 

I applied the latest version of your script i.e. make_site_mutated_protease_v5.py and the latest auxiliary script for splitting the fasta files, along with Zhuofan’s run_v4.1.sh. The script runs but as one might expect, it turned up the following error. I couldn’t find the same error in the log sheet you had prepared.

```bash
ERROR: Assertion `seq_x->ungapped_length() == seq_x->length()` failed.

ERROR:: Exit from: /home/benchmark/rosetta/source/src/core/sequence/Aligner.cc line: 50

Traceback (most recent call last):

 File "../scripts/make_site_mutated_protease_v5.py", line 1032, in <module>

  main(args)

 File "../scripts/make_site_mutated_protease_v5.py", line 976, in main

  ind_type = args.ind_type, pdb_name = args.template_pdb)

 File "../scripts/make_site_mutated_protease_v5.py", line 858, in analyze_mutant_protein

  substitutions, new_subs = compare_sequences(pdb_name, seqrecord.seq, query, ind_type)

 File "../scripts/make_site_mutated_protease_v5.py", line 218, in compare_sequences

  tmpAlign = SWAligner().align(wild_seq_, seq_2_, ss)

RuntimeError:

 

File: /home/benchmark/rosetta/source/src/core/sequence/Aligner.cc:50

[ ERROR ] UtilityExitException

ERROR: Assertion `seq_x->ungapped_length() == seq_x->length()` failed.
```

 

The command I used was : `bash ../scripts/run_v4.1.sh -tp 6m17_chainB_chainE_v3.pdb_relaxed_spikeinterfacev2_11_1.pdb -ml S-protein_GISAID.fasta.txt -ex true -nbh 8.0 -debug true -chk true`



**Reported by**: Changpeng Lu

> 1. For getting mutation results for complex state, I am wondering if we need to do several mutations at a time for each monomer when reference protein we aligned with has mutations in two or three FASTA files.
>
> For example, I am doing with 7bv2, which is a hetrotrimer, and hCoV-19/South_Africa/R05475/2020 | hCoV-19/South_Africa/R05475/2020 reference sequence appears in both Nsp7_GISAID and Nsp8_GISAID. Say if mutation in Nsp7_GISAID for this reference is A123B, while mutation in Nsp8_GISAID for this reference is B234C, should these two mutations be done simultaneously and finally get the mutated Rosetta model that has these two mutations?
>
> ​        Here are replicates in different FASTA files for Nsp7, Nsp8, Nsp12.
>
> ​                
>
> |      | ID 1                             | ID 2                             | in FASTA       |
> | ---- | -------------------------------- | -------------------------------- | -------------- |
> | 1    | hCoV-19/South_Africa/R05475/2020 | hCoV-19/South_Africa/R05475/2020 | Nsp7 and Nsp8  |
> | 2    | hCoV-19/USA/WA-S914/2020         | hCoV-19/USA/WA-S914/2020         | Nsp7 and Nsp8  |
> | 3    | hCoV-19/USA/WA-S1116/2020        | 'hCoV-19/USA/WA-S1116/2020       | Nsp7 and Nsp12 |
> | 4    | hCoV-19/India/S10/2020           | hCoV-19/India/S10/2020           | Nsp8 and Nsp12 |
> | 5    | hCoV-19/USA/NY-NYUMC629/2020     | hCoV-19/USA/NY-NYUMC629/2020     | Nsp8 and Nsp12 |

**Note: All errors could be avoided if using the new script released by Changpeng on July 19th.**



## 2020-07-18 Main Script Updated

**Updated by**: Zhuofan

 I modified the py script so that it will only repack protein residues while hold ligand and RNA motifs fixed using the --only_protein flag. We are now using the same setting to relax the input pdbs.

## 2020-07-19 Main Script New Version Released

**Released by**: Changpeng Lu

**Script Name**: protease_v1.py (later renamed to *make_site_mutated_protein.py*)

Here is a big update for mutation script, named as **“protease_v1.py”.** Now it can handle:

1. Both single multiple FASTA files and multiple FASTA files;

2. In the multiple FASTA case, it could identify the monomer that is corresponded to the specific FASTA files and make mutations for it in complex state;

3. If the reference protein appears in different FASTA files, it will do all mutations read from different FASTA files and make mutations on corresponding monomers, simultaneously;

4. Thanks to Zhuofan, it could handle non-protein part. The basic idea is to make it not moved during repacking.

The command is a little different. Here are several commands with different cases:

1. Multiple FASTA files:

Check out which chains correspond to which FASTA files and use the command like below:



`python make_site_mutated_protease_v5.2.py -t nsp7-nsp8-nsp12_F101A.pdb -m Nsp12_GISAID.fasta.txt, Nsp8_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt -nbh 8.0 -rn test_fast -op -cut "A,B,C" -rn 7bv2_fast`



Here, several new flags are described below:

1) **“-m”** contains all FASTA input files. Use **“,”** **(comma**) to separate them. Make sure not to contain any **space** among them;

2) **“-op”** means avoiding nonprotein residues moving during repacking and minimization;

3) **“-cut”** records the corresponded chain name for each FASTA file. List corresponded chain name and use “,” to separate them. Here, Nsp12 is mapped to Chain A, Nsp8 is mapped to Chain B, Nsp7 is mapped to Chain C. Make sure “-m” and “-cut” flag can map to each other;



2. Single FASTA file:

Command as before:

`  python make_site_mutated_protease_v5.2.py -t nsp7.pdb -m Nsp7_GISAID.fasta.txt -nbh 8.0 -rn test_fast`

**Update about the output:**

​        Since now we need to consider different chains and indices for sequences of chains are not consecutive, I modified the way to show substitutions and mutations as followed:

1) In “mutants.csv”, if “CK7N” in substitutions, it means for chain C, residue 7 is mutated from K in the wild type to N in the reference;

2) In “substitutions.csv”, add a new column “chain” to record the chain name.

Codes have been uploaded to all drive folders, github and attach in this email. Also a small example has been run on 7bv2 and its outputs are attached as well. I am not sure if the updated script could solve Lingjun and Vidur’s problem, but I am wondering if you could run this script for me on your own protein and report errors to me. 

## 2020-07-19 Main Script Debug Fixed

**Reported by**: Lingjun

```bash
ERROR: Error in core::conformation::Conformation::residue(): The sequence position requested was greater than the number of residues in the pose.
ERROR:: Exit from: /scratch/benchmark/W.hojo-2/rosetta.Hojo-2/_commits_/main/source/build/PyRosetta/linux/clang-3.4.2/python-3.5/release/source/src/core/conformation/Conformation.hh line: 508
Traceback (most recent call last):
  File "protease_v1.py", line 1275, in <module>
    main(args)
  File "protease_v1.py", line 1219, in main
    only_protein=args.only_protein)
  File "protease_v1.py", line 1109, in analyze_mutant_protein
    for pm in new_subs]
  File "protease_v1.py", line 1109, in <listcomp>
    for pm in new_subs]
  File "protease_v1.py", line 409, in identify_res_layer
    check_pose = pose.split_by_chain()[main_chain]
RuntimeError: 

File: /scratch/benchmark/W.hojo-2/rosetta.Hojo-2/_commits_/main/source/build/PyRosetta/linux/clang-3.4.2/python-3.5/release/source/src/core/conformation/Conformation.hh:508
[ ERROR ] UtilityExitException
ERROR: Error in core::conformation::Conformation::residue(): The sequence position requested was greater than the number of residues in the pose.
```

The command used below:

`python protease_v1.py -t 6vxx_INPUT.pdb -sym 6vxx.symm -m S-protein_GISAID.fasta.test.txt -nbh 8.0 -rn test_fast`

**Fixed by**: Changpeng Lu

1. For 6vxx, Elliott has done mutations and is relaxing on Amarel. I found the bug from the script is basically on split_by_chain(), that cannot handle symmetrized pose well. So I asked Joey and used an alternative way. Now it worked well for 6vxx case. Once input files are ready, we could get results asap;

## 2020-07-19 Main Script Updated

**Updated by**: Elliott Dolan

The protease_v1.py script is updated to accommodate asymmetric and symmetric membrane proteins. Changpeng is still working on her portion – I have pushed my portion of the code to her repository. When she finishes with her portion, her protease_v1.py upload should contain my modifications. If it doesn’t, I can update it by hand quickly. 

To begin using membrane proteins, add the flag “-memb” or “- - membrane” to declare the protein as a membrane protein. There is nothing intrinsic to a membrane protein pdb file over a soluble protein pdb file, so this argument must be cast.

​        All membrane proteins require a span file, and if you call the membrane protein argument, you must also pass the ‘-mspan’ or ‘—span_file’ argument along with the appropriate span file (it will break if you call -memb without -mspan). If the wrong span file is included, the script will not act properly – it may still run, however. The remainder of the script works well and gives me appropriate .csv outputs.

Additionally, everything for the membrane proteins are uploaded to the CV datafiles folder.

Zhoufan – If you could add something into the run script for membrane proteins, that’d be fantastic. Should be something like this:

Should be callable with `-mem SPAN_FILE`  --- so the -mem flag followed by the span file name

In the case portion of your script:

​       ` -mem) membrane=”$2”;;`

Afterwards:

```bash
If ! [ -z “${membrane}” ]

Then

   Membrane=”-memb -mspan “${membrane}

Fi
```

Finally, just add `${membrane}` onto the python tags for the python script.

Changpeng - The membrane input files have membrane identifying residues, but they don’t seem to throw errors in terms of the pyrosetta alignment code. 

I’ll be around tomorrow to fix things if there are any issues. Hit me up on the slack

## 2020-07-20 Input Data Updated

**Reported by**: Changpeng Lu

> After talking with Zhuofan and Elliott, we found that pdb sequence is not equivalent with the aligned part of wild type sequence in the FASTA file, at least for 6vxx. When we aligned the pdb with the wild type, there are extra mutations except gaps. I am afraid it happens for more proteins that we have done before. In this case, Elliott thought and we agreed the best way to do in the first place is to mutate the pdb with the wild type before relaxing and THEN try to relax and do mutation with those reference proteins in the FASTA. What do you think?

**Updated by**: Elliott Dolan

1. For 6vxx, Elliott has done mutations and is relaxing on Amarel. I found the bug from the script is basically on split_by_chain(), that cannot handle symmetrized pose well. So I asked Joey and used an alternative way. Now it worked well for 6vxx case. Once input files are ready, we could get results asap;

## 2020-07-21 Main Script Debug 

**Reported by**: Changpeng

For 6VSB, the alignment function in Pyrosetta is not a good one. It cannot handle sequences that have many gaps there.

For 6YYT, I got the sense for FASTA files. 

And the script is well running by doubling nsp8 FASTA file in the command below:

`python make_site_mutated_protein.py -t 6yyt.pdb -m Nsp12_GISAID.fasta.txt,Nsp8_GISAID.fasta.txt,Nsp7_GISAID.fasta.txt,Nsp8_GISAID.fasta.txt -nbh 8.0 -op -cut "A,B,C,D" -rn 6yyt_fast`



## 2020-07-22 Main Script Updated

**Updated by**: Changpeng Lu

> I added a new flag called **“-pmm**”, which is equivalent to **“is_pdb_index_match_mutant”**. If the pdb index matches mutant sequence index, the program will align pdb sequence with mutant sequence directly without using pairwise alignment and record the mutation both within and beyond pdb sequence. 

Command below:

`python make_site_mutated_protein.py -t 6VSB_relaxed.pdb -m S-protein_GISAID.fasta.txt,S-protein_GISAID.fasta.txt,S-protein_GISAID.fasta.txt -pmm -nbh 8.0 -rn 6VSB_fast -cut "A,B,C"`

