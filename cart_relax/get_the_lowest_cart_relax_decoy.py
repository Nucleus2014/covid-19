import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-sc', '--scores', type=str, default='score.sc')
args = parser.parse_args()

with open(args.scores, 'r') as scores:
    lines = scores.readlines()

lowest_score = 10000
lowest_decoy = None

for line in lines[2:]:
    info = list(filter(lambda x: x != '', line.split(' ')))
    score = float(info[1])
    decoy = info[-1][:-1]
    if score < lowest_score:
        lowest_score = score
        lowest_decoy = decoy

print(lowest_decoy + '\t' + str(lowest_score))
