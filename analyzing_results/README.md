# commands

## load resuls as list

Concatenate all energy information or point mutations information into a single file as an input of plot scripts down below. Results from different protocols could be concatenated together into one file by setting the same output filename using "-o" flag. Here are general commands that you could use for this script:

**Print out point mutations information into a file**
> python ../../analyzing_results/load_result_as_list.py -i hard-rep_hard-rep/Nsp1_mutants.csv -o Nsp1_results.txt -s 1

**Print out Rosetta total energy scores into a file for make site mutated protein protocol**
> python ../../analyzing_results/load_result_as_list.py -i hard-rep_hard-rep/Nsp1_mutants.csv -o Nsp1_results.txt -s 2

**Cartesian ddG**

> python ../../analyzing_results/load_result_as_list.py -i ddg/Nsp1.ddg -f ddg/Nsp1.fingerprint -o Nsp1_results.txt -s 3
or
> python ../../analyzing_results/load_result_as_list.py -i ddg/Nsp1.ddg -csv hard-rep_hard-rep/Nsp1_mutants.csv -t ddg/Nsp1_cart_INPUT.pdb -o Nsp1_results.txt -s 3

For example, if you would like to concatenate results of all four different protocols for Nsp1 protein. You could follow these steps:

1. Print energy scores for all mutation models of **hard-rep_hard-rep** into the file "Nsp1_results.txt":
> python ../../analyzing_results/load_result_as_list.py -i hard-rep_hard-rep/Nsp1_mutants.csv -o Nsp1_results.txt -s 2
2. Print energy scores for all mutation models of **soft-rep_soft-rep** into the file "Nsp1_results.txt":
> python ../../analyzing_results/load_result_as_list.py -i soft-rep_soft-rep/Nsp1_mutants.csv -o Nsp1_results.txt -s 2
3. Print energy scores for all mutation models of **soft-rep_hard-rep** into the file "Nsp1_results.txt":
> python ../../analyzing_results/load_result_as_list.py -i soft-rep_hard-rep/Nsp1_mutants.csv -o Nsp1_results.txt -s 2
4. Print energy scores for all mutation models of **cartesian ddG** into the file "Nsp1_results.txt":
> python ../../analyzing_results/load_result_as_list.py -i ddg/Nsp1.ddg -f ddg/Nsp1.fingerprint -o Nsp1_results.txt -s 3

## draw plots and histograms
**bar plot**
> python ../../analyzing_results/bar_plot.py -i Nsp1_results.txt -y 4 -t "cart ddG"

**histogram**
> python ../../analyzing_results/histogram.py -i Nsp1_results.txt -data 4 -bins 1 -t "cart ddG"

**scatter plot**
> python ../../analyzing_results/scatter_plot.py -i Nsp1_results.txt -x 2 -y 3 -t "hard-rep+hard-rep vs cart ddG"
