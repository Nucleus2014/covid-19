import os
from statistics import mean

for job in filter(lambda x: os.path.isfile(x + '/' + x + '.fingerprint'), os.listdir()):
    # Get the variant.
    with open(job + '/' + job + '.fingerprint', 'r') as p_fingerprint:
        variant_name = p_fingerprint.read()[:-1]
    # If there are point mutations in the PDB model.
    if os.path.isfile(job + '/' + job + '.ddg'):
        # Read scores from the ddG file.
        with open(job + '/' + job + '.ddg', 'r') as p_ddg:
            ddg_lines = p_ddg.readlines()
        # Obtain the wild type score.
        wt_decoy_scores = list()
        for i in range(3):
            wt_decoy_scores.append(float(list(filter(lambda x: x != '', ddg_lines[i].split(' ')))[3]))
        wt_score = mean(wt_decoy_scores)
        # Obtain the variant score.
        variant_decoy_scores = list()
        for i in range(3, 6):
            variant_decoy_scores.append(float(list(filter(lambda x: x != '', ddg_lines[i].split(' ')))[3]))
        variant_score = mean(variant_decoy_scores)
        # Calculate the ddG value.
        variant_ddg = str(round(variant_score - wt_score, 2))
        # Convert the ddg file to a ddg2 file.
        with open(job + '/' + job + '.csv', 'w+') as p_csv:
            p_csv.write(variant_name + ',' + variant_ddg + '\n')
    # If no point mutations present in the PDB model.
    else:
        with open(job + '/' + job + '.csv', 'w+') as p_csv:
            p_csv.write(variant_name + ',NA\n')
