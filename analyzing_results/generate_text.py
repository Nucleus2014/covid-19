import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-idx', '--line_index', type=int, required=True)
parser.add_argument('-thld', '--threshold', type=float, default=10.0)
parser.add_argument('-o', '--output_prefix', type=str)
args = parser.parse_args()

with open(args.input_file, 'r') as p_input:
    lines = p_input.readlines()

# In the "lines" list, line number starts at 0 instead of 1.
ddg_values = eval(lines[args.line_index - 1])

non_zero_total_num = 0
more_stablized_total_num = 0
more_stablized_total_score = 0
less_stablized_total_num = 0
less_stablized_total_score = 0
outlier_total_num = 0

for ddg in ddg_values:
    if ddg > 0:
        non_zero_total_num += 1
        less_stablized_total_num += 1
        less_stablized_total_score += ddg
        if ddg > args.threshold:
            outlier_total_num += 1
    elif ddg < 0:
        non_zero_total_num += 1
        more_stablized_total_num += 1
        more_stablized_total_score += ddg

more_stablized_percent = round(more_stablized_total_num / non_zero_total_num * 100, 2)
less_stablized_percent = round(less_stablized_total_num / non_zero_total_num * 100, 2)
outlier_percent = round(outlier_total_num / non_zero_total_num * 100, 2)

more_stablized_average_score = round(-more_stablized_total_score / more_stablized_total_num, 2)
less_stablized_average_score = round(less_stablized_total_score / less_stablized_total_num, 2)

if not args.output_prefix:
    args.output_prefix = args.input_file[:args.input_file.find('_')]

with open(args.output_prefix + '.txt', 'w+') as pf:
    string = 'The vast majority of the USVs ({}%) were estimated to be slightly less stable than the reference sequence (average change in apparent free energy of stabilization ~+{} REUs). In less than {}% cases, the estimated change in apparent free energy of stabilization exceeded +{} REUs, which we attribute to incomplete sampling of alternative backbone atom positions following introduction of the amino acid change(s) in the homology model. A minority of the USVs ({}%) were estimated to be more stable than the reference sequence (average change in apparent free energy of stabilization ~-{} REUs).\n'.\
        format(str(less_stablized_percent), str(less_stablized_average_score), \
        str(outlier_percent), str(args.threshold), str(more_stablized_percent), \
        str(more_stablized_average_score))
    pf.write(string)
