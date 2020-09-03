import argparse, sys
from Bio import SeqUtils
from pyrosetta import *


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True, help='*.ddg or *mutant.csv')
    parser.add_argument('-csv', '--csv_file', type=str, help='for *.ddg input file')
    parser.add_argument('-t', '--template_pdb', type=str, help='for *.ddg input file')
    parser.add_argument('-f', '--fingerprint_file', type=str, help='for *.ddg input file')
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-s', '--step', type=int, choices=[1,2,3], required=True)
    return parser.parse_args()

def generate_scores_dict_from_ddg(ddg):
    variants = dict()
    with open(ddg, 'r') as p_ddg:
        for line in p_ddg:
            data = list(filter(lambda x: x != '', line.split(' ')))
            name = data[0]
            ddg = data[1][:-1]
            variants[name] = ddg
    return variants

def generate_list_from_csv_pdb(input_file, csv, pdb):
    info = '['
    init()
    pose = pose_from_pdb(pdb)
    with open(csv, 'r') as p2csv:
        csv_lines = p2csv.readlines()
    for idx, item in enumerate(lines[0].split(',')):
        if item == 'substitutions':
            break
    # Calculate ddG
    ddg_variants_dict = generate_scores_dict_from_ddg(input_file)
    for csv_line in csv_lines[1:]:
        mutations = csv_line.split(',')[idx] # 'AI559_V;AP585_S'
        if mutations == '':
            info += '0,'
        else:
            key = 'MUT'
            for mutation in mutations.split(';'):
                chain_id = mutation[0]
                # native_aa = mutation[1]
                mutated_res = SeqUtils.IUPACData.protein_letters_1to3[mutation[-1]].upper()
                pdb_index = int(mutation[2:-2])
                pose_index = pose.pdb_info().pdb2pose(chain_id, pdb_index)
                key += '_' + pose_index + mutated_res
            ddg = ddg_variants_dict[key]
            info += str(round(ddg, 2)) + ','
    info = info[:-1] + ']\n'
    return info

def generate_list_from_fingerprint(input_file, fingerprint):
    info = '['
    # Calculate ddG
    ddg_variants_dict = generate_scores_dict_from_ddg(input_file)
    with open(fingerprint, 'r') as p2fingerprint:
        for line in p2fingerprint:
            if line.startswith('WT'):
                info += '0,'
            else:
                key = 'MUT'
                mutations = line[:-1].split(',')
                mut_order = {}
                for mutation in mutations:
                    mutantion_info = mutation.split(' ')
                    pose_index = mutantion_info[1]
                    mutated_res = SeqUtils.IUPACData.protein_letters_1to3[mutantion_info[2]].upper()
                    mut_order[pose_index] = mutated_res
                keys = [int(x) for x in mut_order.keys()]
                for keys in sorted(keys):
                    key += '_' + str(keys) + mut_order[str(keys)] 
                info += ddg_variants_dict[key] + ','
    info = info[:-1] + ']\n'
    return info


args = parse_arguments()
if args.step == 1:
    with open(args.input_file, 'r') as p_input:
        lines = p_input.readlines()
    info = '['
    for idx, item in enumerate(lines[0].split(',')):
        if item == 'substitutions':
            break
    for line in lines[1:]:
        mut = line.split(',')[idx][1:]
        if mut == '':
            mut = 'WT'
        info += '"' + mut + '",'
    info = info[:-1] + ']\n'
elif args.step == 2:
    with open(args.input_file, 'r') as p_input:
        lines = p_input.readlines()
    info = '['
    for idx, item in enumerate(lines[0].split(',')):
        if item == 'energy_change':
            break
    for line in lines[1:]:
        ddg = line.split(',')[idx]
        if ddg == 'NA':
            ddg = '0.00'
        info += str(round(float(ddg), 2)) + ','
    info = info[:-1] + ']\n'
elif args.step == 3:
    if args.csv_file:
        info = generate_list_from_csv_pdb(args.input_file, args.csv_file, args.template_pdb)
    elif args.fingerprint_file:
        info = generate_list_from_fingerprint(args.input_file, args.fingerprint_file)

with open(args.output_file, 'a+') as report:
    report.write(info)
