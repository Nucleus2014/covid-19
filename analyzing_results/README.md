python ..\..\analyze_results\bar_plot.py -i .\Nsp1_results.txt -y 4 -t "cart ddG"
python ..\..\analyze_results\histogram.py -i .\Nsp1_results.txt -data 4 -bins 1 -t "cart ddG"
python ..\..\analyze_results\scatter_plot.py -i Nsp1_results.txt -t "hh vs cart" -x 2 -y 3