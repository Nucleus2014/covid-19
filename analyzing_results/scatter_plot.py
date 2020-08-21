import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-t', '--title', type=str)
parser.add_argument('-x', '--x_axis', type=int, nargs='*', required=True)
parser.add_argument('-y', '--y_axis', type=int, nargs='*', required=True)
parser.add_argument('-c', '--colors', type=str, nargs='*')
parser.add_argument('-m', '--markers', type=str, nargs='*')
args = parser.parse_args()


with open(args.input_file, 'r') as p_input:
    lines = p_input.readlines()

for i in range(len(args.x_axis)):
    x_values = eval(lines[args.x_axis[i] - 1]) # In the "lines" list, line number starts at 0 instead of 1.
    y_values = eval(lines[args.y_axis[i] - 1])
    if args.colors:
        color = args.colors[i]
    else:
        color = 'b'
    if args.markers:
        marker = args.markers[i]
    else:
        marker = None
    plt.scatter(x_values, y_values, s=5, c=color, marker=marker)

if args.title:
    plt.title(args.title)

# Always remember that points beyond the scope will not be included
axes = plt.gca()
axes.set_xlim([-6, 10])
axes.set_ylim([-6, 10])

fig = plt.gcf()
fig.set_size_inches(20, 20)

plt.show()
# fig.savefig('test2png.png', dpi=100)
