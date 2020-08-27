import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-data', '--line_index', type=int, required=True)
parser.add_argument('-bins', '--bins', type=float, default=1)
parser.add_argument('-t', '--title', type=str)
parser.add_argument('-o', '--output_prefix', type=str)
args = parser.parse_args()

with open(args.input_file, 'r') as p_input:
    lines = p_input.readlines()

# In the "lines" list, line number starts at 0 instead of 1.
data = eval(lines[args.line_index - 1])

plt.hist(data, bins=args.bins)

plt.xlabel('ddG')
plt.ylabel('Numbers of variants')

plt.title(args.title)

fig = plt.gcf()
# fig.set_size_inches(10, 8)

plt.show()

if args.output_prefix:
    fig.savefig(args.output_prefix + '.png', dpi=100)
