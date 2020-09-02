while (( $# > 1 ))
do
    case $1 in
        -p) protein_name="$2";;
        *) break;
    esac; shift 2
done

python ../../analyzing_results/bar_plot.py -i ${protein_name}_results.txt -y 2 -t "${protein_name} hard-rep hard-rep" -o "bar_plots_${protein_name}_hard-rep_hard-rep"
python ../../analyzing_results/bar_plot.py -i ${protein_name}_results.txt -y 3 -t "${protein_name} hard-rep soft-rep" -o "bar_plots_${protein_name}_hard-rep_soft-rep"
python ../../analyzing_results/bar_plot.py -i ${protein_name}_results.txt -y 4 -t "${protein_name} soft-rep soft-rep" -o "bar_plots_${protein_name}_soft-rep_soft-rep"
# python ../../analyzing_results/bar_plot.py -i ${protein_name}_results.txt -y 5 -t ${protein_name}" Cart ddG" -o "bar_plots_"${protein_name}"_Cart_ddG"
