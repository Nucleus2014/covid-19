import os
from statistics import mean

for ddg in filter(lambda file: file.endswith('.ddg'), os.listdir()):
    with open(ddg, 'r') as p_ddg:
        lines = p_ddg.readlines()
    # Obtain the wild type score
    wt_decoy_scores = list()
    for i in range(3):
        wt_decoy_scores.append(float(list(filter(lambda x: x != '', lines[i].split(' ')))[3]))
    wt_score = mean(wt_decoy_scores)
    # Convert the ddg file to a ddg2 file
    with open(ddg + '2', 'w+') as p_ddg2:
        for i in range(3, len(lines), 3):
            # Obtain the mutations
            variant_name = list(filter(lambda x: x != '', lines[i].split(' ')))[2][:-1]
            # Obtain the variant score
            variant_decoy_scores = list()
            for j in range(3):
                variant_decoy_scores.append(float(list(filter(lambda x: x != '', lines[i + j].split(' ')))[3]))
            variant_score = mean(variant_decoy_scores)
            # Calculate ddG value
            variant_ddg = str(round(variant_score - wt_score, 2))
            p_ddg2.write(variant_name + ' ' + variant_ddg + '\n')
