# This script is to check whether reference proteins are the same or not for different GISAID files
# Editor: Changpeng Lu
# Date: 2020-07-16

import os
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', type=str, help="list all fasta files that needs to be compared")
args = parser.parse_args()

def get_id(fasta_name):
    fp = open(fasta_name,"r")
    labels = []
    for line in fp:
        if line[0] == ">":
            labels.append(line.strip()[1:].split("|")[3:5])
    labels.pop(0)
    fp.close()
    return labels

list_fasta_names = args.input.strip().split(",")
print("------------------------------")
print("Considered FASTA files are:")
print(list_fasta_names)
ids = []
for name in list_fasta_names:
    id_1_2 = get_id(name)
    ids.append(id_1_2)

# provide mapping between several fasta files
preset = []
indices = []
for i in range(len(ids)):
    tmp = list(range(len(ids)))
    tmp.pop(i)
    for j in tmp:
        if (i,j) not in preset:
            print(preset)
            preset.append((i,j))
            preset.append((j,i))
            indices.append((i,j))

replicates = []
for mapping in indices:
    for que in ids[mapping[0]]:
        try:
            tmp_ind = ids[mapping[1]].index(que)
            replicates.append((que,list_fasta_names[mapping[0]],list_fasta_names[mapping[1]]))
        except ValueError:
            pass
for j in replicates:
    print("\n-----------------------------\n")
    print("ID1 and ID2 for the replicates protein are:")
    print(j[0])
    print("It is in both files:{0},{1}".format(j[1],j[2]))
    print("\n")
