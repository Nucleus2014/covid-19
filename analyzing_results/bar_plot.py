import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-y', '--y_line_index', type=int, required=True)
parser.add_argument('-x', '--x_line_index', type=int, help='Optional. Passing this flag will show the point mutations at the x axis.')
parser.add_argument('-t', '--title', type=str)
parser.add_argument('-o', '--output_prefix', type=str)
args = parser.parse_args()

with open(args.input_file, 'r') as p_input:
    lines = p_input.readlines()

# In the "lines" list, line number starts at 0 instead of 1.
y_values = eval(lines[args.y_line_index - 1])

x = list(range(len(y_values)))
plt.bar(x, y_values, width=0.5, color='blue')

if args.x_line_index:
    x_ticks = eval(lines[args.x_line_index - 1])
    plt.xticks(x, x_ticks)
else:
    plt.xticks(x, ['' for i in x])

# plt.yticks(list(range(-7, 11)), list(range(-7, 11)))

plt.xlabel('Variants')
plt.ylabel('ddG')

plt.title(args.title)\

fig = plt.gcf()
# fig.set_size_inches(10, 8)

plt.show()

if args.output_prefix:
    fig.savefig(args.output_prefix + '.png', dpi=100)
