## commands

###### load resuls as list
**Point mutations information**
> python ../../analyzing_results/load_result_as_list.py -i hard-rep_hard-rep/Nsp1_mutants.csv -o Nsp1_results.txt -s 1
**Make site mutated protein protocol**
> python ../../analyzing_results/load_result_as_list.py -i hard-rep_hard-rep/Nsp1_mutants.csv -o Nsp1_results.txt -s 2
**Cartesian ddG**
If you use *.fasta files to generate *.mut files and *.fingerprint files:
> python ../../analyzing_results/load_result_as_list.py -i ddg/Nsp1.ddg -f ddg/Nsp1.fingerprint -o Nsp1_results.txt -s 3
If you use *_mutants.csv file to generate *.mut files:
> python ../../analyzing_results/load_result_as_list.py -i ddg/Nsp1.ddg -csv hard-rep_hard-rep/Nsp1_mutants.csv -o Nsp1_results.txt -s 3
If you use *_mutants.csv file (probably with some previous output *.ddg files) to generate *.mut files, but the results in the output *.ddg file is not in the correct order:
> python ../../analyzing_results/load_result_as_list.py -i ddg/Nsp1.ddg -csv hard-rep_hard-rep/Nsp1_mutants.csv -t ddg/Nsp1_cart_INPUT.pdb -o Nsp1_results.txt -s 3

###### draw plots and histograms
**bar plot**
> python ../../analyzing_results/bar_plot.py -i Nsp1_results.txt -y 4 -t "cart ddG"
**histogram**
> python ../../analyzing_results/histogram.py -i Nsp1_results.txt -data 4 -bins 1 -t "cart ddG"
**scatter plot**
> python ../../analyzing_results/scatter_plot.py -i Nsp1_results.txt -x 2 -y 3 -t "hard-rep+hard-rep vs cart ddG"
