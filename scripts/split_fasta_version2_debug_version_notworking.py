# This script is to split fasta file into n splits
# Editor: Changpeng Lu
# Date: 2020/06/30

# Usage:
'''
	python split_fasta.py -i <fasta-file-name> -n <number of splits, default to 12>
'''

# Initial startup
import os
import warning
import argparse
import math
import pyrosetta as pr
from pyrosetta.rosetta.core.sequence import Aligner, Sequence, read_fasta_file
from pyrosetta.io import pose_from_sequence
#from Bio import Align


# All flags
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, help='the filename of fasta file')
parser.add_argument('-n', '--splits', type=int, default=12,
                    help='number of splits you want to divide fasta file. Suggestion is to keep the same with the number of nodes that is offered on the server')
parser.add_argument('-t', '--template_pdb',type=str)
args = parser.parse_args()

pr.init()
n = args.splits
filename = args.input

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
split_size = math.floor(len(seqs) / n)
query_value = seqs.pop(query_name)
keys_list.pop(0)
print("-----------------------------------------")
print("The number of reference sequences are: {}".format(len(keys_list)))
print("Split into {} pieces".format(n))
print("Each piece is about to have {} number of reference sequences.".format(split_size))

wild_pose = pr.pose_from_pdb(args.template_pdb)
wild_seq = wild_pose.chain_sequence(1)
if len(wild_seq) != len(query_value): # pdb not match length of query in fasta file
	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		warnings.warn("RUNTIME WARNING",RuntimeWarning)
	#print("WARNING")
		print("QUERY not match with pdb! Split files are all adjusted. However, query sequence in the fasta file still not the right one. Please run this code again with -n 1")
	seqs_set = read_fasta_file(args.input)
	refs = []
	n = 0
	for id in keys_list:
		if n == 0:
			query = seqs[id]
		refs.append(Sequence(keys_list[id],id))
		n += 1
	wild_seq_ = Sequence(wild_seq,"pdb")
	refs_out = []
	for s in refs:
		ss = SimpleScoringScheme()
		tmpAlign = SWAligner().align(wild_seq_,s,ss)
		ref_out = tmpAlign.to_string().strip().split("\n")[2].strip().split()
		refs_out.append(ref_out[2]) # save reference sequence
	#aligner = Align.PairwiseAligner()
	#aligner.mode = 'global'
	#aligner.match_score = 2
	#aligner.mismatch_score = -1
	#alignments = aligner.align(wild_seq, query_value)
	aligner = Aligner()
	#seq = Sequence()
	#seq.append_char(wild_seq)
	#qseq = Sequence()
	#qseq.append_char(query_value)
	#aa = aligner.align(seq, qseq)
	#print(aa)
	seqs_set = read_fasta_file(args.input)
	aa = aligner.align(seqs_set[0],seqs_set[1])
	#query_value = wild_seq
else:
	print("QUERY match with wild_seq")
wrongs = []
for j in range(len(keys_list)):
	if len(seqs[keys_list[j]]) != len(query_value):
		print("WARNING:\n")
		print("Length of Sequence not matched! Reference sequences are being adjusted based on query sequence.")
		print(seqs[keys_list[j]])
		print(len(seqs[keys_list[j]]))
		wrongs.append(keys_list[j]) # save keys of wrong seqs	
if len(wrongs) != 0:
	print("The number of sequences that should be modified is:{}".format(len(wrongs)))
	for each in wrongs:
		tmp = seqs[each]
		start = 0
		while query_value[start] != tmp[start]:
			start += 1
		new_tmp = tmp[start:]
		if len(new_tmp) != len(query_value):
			new_tmp = new_tmp[:len(query_value)]
		seqs[each] = new_tmp
	print(new_tmp)
	print("-----------------------")
	print("Protein sequence is adjusted to the sequence in PDB file:\n{}".format(query_value))
	#print(aa)
for i in range(n):
	if n == 1:
		gp = open(args.input + ".modified.txt", "w")
	else:
		gp = open(args.input + "." + str(i+1) + ".txt", "w")
	gp.write(query_name + "\n")
	seq_split_batch = math.floor(len(query_value) / 60)
	for batch in range(seq_split_batch):
		tmp = query_value[batch * 60 : (batch + 1) * 60]
		gp.write(str(tmp) + "\n")
	gp.write(query_value[(batch + 1) * 60 :] + "\n")	
	#for j in query_value:
	#	gp.write(j + "\n")
	if n == 1:
		for j in range(len(keys_list)):
			gp.write(keys_list[j] + "\n")
			seq_split_batch = math.floor(len(seqs[keys_list[j]]) / 60)
			for batch in range(seq_split_batch):
				tmp = seqs[keys_list[j]][batch * 60 : (batch + 1) * 60]
				gp.write(str(tmp) + "\n")
			gp.write(seqs[keys_list[j]][(batch + 1) * 60 :] + "\n")
		break
	if i < n -1:
		for j in range(i * split_size, (i+1) * split_size):
			gp.write(keys_list[j] + "\n")
			seq_split_batch = math.floor(len(seqs[keys_list[j]]) / 60)
			for batch in range(seq_split_batch):
				tmp = seqs[keys_list[j]][batch * 60 : (batch + 1) * 60]
				gp.write(str(tmp) + "\n")
			gp.write(seqs[keys_list[j]][(batch + 1) * 60 :] + "\n")
			#for s in seqs[keys_list[j]]:
				#gp.write(s + "\n")
	else:
		if len(seqs) % n != 0:
			remain_keys = keys_list[i * split_size:]
			for j in range(0,len(remain_keys)):
				gp.write(remain_keys[j] + "\n")
				seq_split_batch = math.floor(len(seqs[keys_list[j]]) / 60)
				for batch in range(seq_split_batch):
					gp.write(seqs[keys_list[j]][batch * 60 : (batch + 1) * 60])
					gp.write("\n")
				gp.write(seqs[keys_list[j]][(batch + 1) * 60 :])
				#for s in seqs[remain_keys[j]]:
					#gp.write(s + "\n")
	gp.close()

