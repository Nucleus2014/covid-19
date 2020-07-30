# This script is to split fasta file into n splits
# Editor: Changpeng Lu
# Date: 2020/06/30

# Usage:
'''
	python split_fasta.py -i <fasta-file-name> -n <number of splits, default to 12> --check
'''

# Initial startup
import os
import warnings
import argparse
import math
import pyrosetta as pr
from pyrosetta.rosetta.core.sequence import SWAligner, Sequence, read_fasta_file,\
SimpleScoringScheme
from pyrosetta.io import pose_from_sequence
import numpy as np


# All flags
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, help='the filename of fasta file')
parser.add_argument('-n', '--splits', type=int, default=12,
                    help='number of splits you want to divide fasta file. Suggestion is to keep the same with the number of nodes that is offered on the server')
parser.add_argument('-t', '--template_pdb',type=str)
parser.add_argument('--check', action='store_true',help="do pairwise sequence alignment and adjust both separate and original fasta file")
args = parser.parse_args()

pr.init()
n_splits = args.splits
filename = args.input

# functions
def find(s, ch): # find all characters and return their indices
	return [i for i, ltr in enumerate(s) if ltr == ch] 

seqs = {}
fp = open(filename,'r')
c = 0
count = 0
keys_list = []
for line in fp:
	if line[0] == ">":
		if count == 0:
			query_name = line.strip()
		c += 1
		if c > 1:
			seqs[tag] = seq
			c = 1
		seq = ""
		tag = line.strip()
		keys_list.append(tag)
		count += 1
	else:
		#seq.append(line.strip())
		seq += line.strip()
seqs[tag] = seq # save the last one
fp.close()
query_value = seqs.pop(query_name)
keys_list.pop(0)
split_size = math.ceil(len(seqs) / n_splits)
refs_out = [seqs[kk] for kk in keys_list]
print("-----------------------------------------")
print("The number of reference sequences are: {}".format(len(keys_list)))
print("Split into {} pieces".format(n_splits))
print("Each piece is about to have {} number of reference sequences.".format(split_size))

wild_pose = pr.pose_from_pdb(args.template_pdb)
wild_seq = wild_pose.chain_sequence(1)
print("sequence in pdb is:{}".format(wild_seq))

if len(wild_seq) != len(query_value): # pdb not match length of query in fasta file
	warnings.simplefilter("always")
	warnings.warn("WARNING:QUERY not match with pdb! Adjust automatically by using --check flag",RuntimeWarning)

if np.sum([len(x) != len(wild_seq) for x in seqs.values()]) != 0:
	warnings.simplefilter("always")
	warnings.warn("The length of mutants' sequences not match to the length of the protein sequence in pdb file!", RuntimeWarning)
if args.check == True:
	#seqs_set = read_fasta_file(filename)
	refs = []
	n = 0
	if "|" in query_name:
		tmp_query_name = query_name.split("|")
		query_name = tmp_query_name[0]
	refs.append(Sequence(query_value, query_name))
	for id in keys_list:
		refs.append(Sequence(seqs[id],id)) 
		n += 1
	wild_seq_ = Sequence(wild_seq,"pdb")
	refs_out = []
	n = 0
	for s in refs:
		ss = SimpleScoringScheme()
		tmpAlign = SWAligner().align(wild_seq_,s,ss)
		start_ind = int(tmpAlign.to_string().strip().split("\n ")[1].split()[1])
		if start_ind != 1:
			hyps = wild_seq[0: start_ind - 1]
		else:
			hyps = ""
		if n == 0:
			que_out = tmpAlign.to_string().strip().split("\n ")[2].strip().split(" ")[-1].replace("\n", "")
			que_out = hyps + que_out
		else:
			ref_out = tmpAlign.to_string().strip().split("\n>")[1].strip().split(" ")[-1].replace("\n", "")
			ref_out = hyps + ref_out
			refs_out.append(ref_out)		
		n += 1
	print("After adjusting, the number of residues for sequences in fasta match with the number of residues, which is:{}".format(len(wild_seq)))
	print("NOTE: This number of residues may not be the maximum number of the pose, since starting index in pose may be larger than 1.")
	if len(wild_seq) == len(que_out.strip()):
		print("length of pdb sequence is equal to length of query sequence!")		
else:
	print("QUERY match with wild_seq")

for ind in range(n_splits):
	if args.check == False:
		que_out = query_value
	gp = open(args.input + "." + str(ind+1) + ".txt", "w")
	gp.write(query_name + "\n")
	seq_split_batch = math.floor(len(que_out) / 60)
	for batch in range(seq_split_batch):
		tmp = que_out[batch * 60 : (batch + 1) * 60]
		gp.write(str(tmp) + "\n")
	if seq_split_batch != 0:
		gp.write(que_out[(batch + 1) * 60 :] + "\n")
	else:
		gp.write(que_out + "\n")	
	
	if n_splits == 1:
		for j in range(len(keys_list)):
			gp.write(keys_list[j] + "\n")
			seq_split_batch = math.floor(len(refs_out[j]) / 60)
			for batch in range(seq_split_batch):
				tmp = refs_out[j][batch * 60 : (batch + 1) * 60]
				gp.write(str(tmp) + "\n")
			if seq_split_batch != 0:
				gp.write(refs_out[j][(batch + 1) * 60 :] + "\n")
			else:
				gp.write(refs_out[j] + "\n")
		break
	if ind < n_splits -1:
		for k in range(ind * split_size, (ind+1) * split_size):
			gp.write(keys_list[k] + "\n")
			seq_split_batch = math.floor(len(refs_out[k]) / 60)
			for batch in range(seq_split_batch):
				tmp = refs_out[k][batch * 60 : (batch + 1) * 60]
				gp.write(str(tmp) + "\n")
			if seq_split_batch != 0:
				gp.write(refs_out[k][(batch + 1) * 60 :] + "\n")
			else:
				gp.write(refs_out[k] + "\n")
	else:
		remain_keys = keys_list[ind * split_size:]
		remain_refs = refs_out[ind * split_size:]
		for q in range(0,len(remain_keys)):
			gp.write(remain_keys[q] + "\n")
			seq_split_batch = math.floor(len(remain_refs[q]) / 60)
			for batch in range(seq_split_batch):
				gp.write(remain_refs[q][batch * 60 : (batch + 1) * 60])
				gp.write("\n")
			if seq_split_batch != 0:
				gp.write(remain_refs[q][(batch + 1) * 60 :] + "\n")
			else:
				gp.write(remain_refs[q] + "\n")
	gp.close()

if args.check == True:
	cp = open(args.input + ".modified.txt", "w")
	cp.write(query_name + "\n")
	seq_split_batch = math.floor(len(que_out) / 60)
	for batch in range(seq_split_batch):
		tmp = que_out[batch * 60 : (batch + 1) * 60]
		cp.write(str(tmp) + "\n")
	if seq_split_batch != 0:
		cp.write(que_out[(batch + 1) * 60 :] + "\n")
	else:
		cp.write(que_out + "\n")	

	for p in range(len(keys_list)):
		cp.write(keys_list[p] + "\n")
		seq_split_batch = math.floor(len(refs_out[p]) / 60)
		for batch in range(seq_split_batch):
			tmp = refs_out[p][batch * 60 : (batch + 1) * 60]
			cp.write(str(tmp) + "\n")
		if seq_split_batch != 0:
			cp.write(refs_out[p][(batch + 1) * 60 :] + "\n")
		else:
			cp.write(refs_out[p] + "\n")
	cp.close()
