while (( $# > 1 ))
do
    case $1 in
        -p) protein_name="$2";;
        *) break;
    esac; shift 2
done

python ../../analyzing_results/histogram.py -i ${protein_name}_results.txt -data 2 -bins 10 -t "${protein_name} hard-rep hard-rep" -o "histograms_${protein_name}_hard-rep_hard-rep"
python ../../analyzing_results/histogram.py -i ${protein_name}_results.txt -data 3 -bins 10 -t "${protein_name} soft-rep hard-rep" -o "histograms_${protein_name}_soft-rep_hard-rep"
python ../../analyzing_results/histogram.py -i ${protein_name}_results.txt -data 4 -bins 10 -t "${protein_name} soft-rep soft-rep" -o "histograms_${protein_name}_soft-rep_soft-rep"
python ../../analyzing_results/histogram.py -i ${protein_name}_results.txt -data 5 -bins 10 -t "${protein_name} Cart ddG" -o "histograms_${protein_name}_Cart_ddG"
python ../../analyzing_results/histogram.py -i ${protein_name}_results.txt -data 6 -bins 10 -t "${protein_name} FastRelax" -o "histograms_${protein_name}_FastRelax"
