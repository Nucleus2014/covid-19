while (( $# > 1 ))
do
    case $1 in
        -p) protein_name="$2";;
        *) break;
    esac; shift 2
done

python ../../analyzing_results/scatter_plot.py -i ${protein_name}_results.txt -x 5 -y 2 -t "${protein_name} hard-rep hard-rep vs Cart ddG" -o "scatter_plots_${protein_name}_hard-rep_hard-rep_vs_Cart_ddG"
python ../../analyzing_results/scatter_plot.py -i ${protein_name}_results.txt -x 5 -y 3 -t "${protein_name} soft-rep hard-rep vs Cart ddG" -o "scatter_plots_${protein_name}_soft-rep_hard-rep_vs_Cart_ddG"
python ../../analyzing_results/scatter_plot.py -i ${protein_name}_results.txt -x 5 -y 4 -t "${protein_name} soft-rep soft-rep vs Cart ddG" -o "scatter_plots_${protein_name}_soft-rep_soft-rep_vs_Cart_ddG"
python ../../analyzing_results/scatter_plot.py -i ${protein_name}_results.txt -x 5 -y 6 -t "${protein_name} FastRelax vs Cart ddG" -o "scatter_plots_${protein_name}_FastRelax_vs_Cart_ddG"
