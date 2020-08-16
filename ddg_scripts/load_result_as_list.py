import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-f', '--fingerprint_file', type=str)
parser.add_argument('-o', '--output_file', type=str, required=True)
parser.add_argument('-s', '--step', type=int, required=True)
args = parser.parse_args()

with open(args.input_file, 'r') as p_input:
    lines = p_input.readlines()
info = '['
if args.step == 1:
    for idx, item in enumerate(lines[0].split(',')):
        if item == 'substitutions':
            break
    for line in lines[1:]:
        mut = line.split(',')[idx][1:]
        if mut == '':
            mut = 'WT'
        info += '"' + mut + '",'
elif args.step == 2:
    for idx, item in enumerate(lines[0].split(',')):
        if item == 'energy_change':
            break
    for line in lines[1:]:
        ddg = line.split(',')[idx]
        if ddg == 'NA':
            ddg = '0.00'
        info += str(round(float(ddg), 2)) + ','
elif args.step == 3:
    # Get the wild type score
    wt_scores = list()
    for i in range(3):
        wt_scores.append(float(list(filter(lambda x : x != '', lines[i].split(' ')))[3]))
    wt_score = min(wt_scores)
    # Obtain scores of each variant
    variants_scores = list()
    for i in range(3, len(lines), 3):
        scores = list()
        for j in range(3):
            scores.append(float(list(filter(lambda x : x != '', lines[i + j].split(' ')))[3]))
        variants_scores.append(min(scores))
    # Calculate ddG
    diff = 0
    with open(args.fingerprint_file, 'r') as p_fingerprint:
        for idx, line in enumerate(p_fingerprint):
            if line == 'WT,\n':
                info += '0,'
                diff += 1
            else:
                info += str(round(variants_scores[idx - diff] - wt_score, 2)) + ','
info = info[:-1] + ']\n'
with open(args.output_file, 'a+') as report:
    report.write(info)
