lowest_score = 10000
lowest_decoy = None

for i in range(1, 101):
    with open(str(i) + '/score.sc', 'r') as scores:
        line = scores.readlines()[-1]
    info = list(filter(lambda x: x != '', line.split(' ')))
    score = float(info[1])
    if score < lowest_score:
        lowest_score = score
        lowest_decoy = i

print(str(lowest_decoy) + '\t' + str(lowest_score))
