import argparse
from Bio import SeqUtils
# from pyrosetta import *


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
    with open(ddg, 'r') as p2ddg:
        lines = p2ddg.readlines()
    # Obtain the wild type score
    wt_decoy_scores = list()
    for i in range(3):
        wt_decoy_scores.append(float(list(filter(lambda x: x != '', lines[i].split(' ')))[3]))
    wt_score = min(wt_decoy_scores)
    # Obtain scores of each variant
    variants = dict()
    for i in range(3, len(lines), 3):
        variant_name = list(filter(lambda x: x != '', lines[i].split(' ')))[2][:-1]
        decoy_scores = list()
        for j in range(3):
            decoy_scores.append(float(list(filter(lambda x: x != '', lines[i + j].split(' ')))[3]))
        ddg = min(decoy_scores) - wt_score
        variants[variant_name] = ddg
    return variants

def generate_scores_list_from_ddg(ddg):
    with open(ddg, 'r') as p2ddg:
        lines = p2ddg.readlines()
    # Obtain the wild type score
    wt_decoy_scores = list()
    for i in range(3):
        wt_decoy_scores.append(float(list(filter(lambda x: x != '', lines[i].split(' ')))[3]))
    wt_score = min(wt_decoy_scores)
    # Obtain scores of each variant
    variants = list()
    for i in range(3, len(lines), 3):
        decoy_scores = list()
        for j in range(3):
            decoy_scores.append(float(list(filter(lambda x: x != '', lines[i + j].split(' ')))[3]))
        ddg = min(decoy_scores) - wt_score
        variants.append(ddg)
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

def generate_list_from_csv(input_file, csv):
    info = '['
    with open(csv, 'r') as p2csv:
        csv_lines = p2csv.readlines()
    for idx, item in enumerate(csv_lines[0].split(',')):
        if item == 'substitutions':
            break
    # Calculate ddG
    ddg_variants_list = generate_scores_list_from_ddg(input_file)
    diff = 0
    for i, csv_line in enumerate(csv_lines[1:]):
        mutations = csv_line.split(',')[idx] # 'AI559_V;AP585_S'
        if mutations == '':
            info += '0,'
            diff += 1
        else:
            info += str(round(ddg_variants_list[i - diff], 2)) + ','
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
                for mutation in mutations:
                    mutantion_info = mutation.split(' ')
                    pose_index = mutantion_info[1]
                    mutated_res = SeqUtils.IUPACData.protein_letters_1to3[mutantion_info[2]].upper()
                    key += '_' + pose_index + mutated_res
                ddg = ddg_variants_dict[key]
                info += str(round(ddg, 2)) + ','
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
        if args.pdb:
            info = generate_list_from_csv_pdb(args.input_file, args.csv_file, args.template_pdb)
        else:
            info = generate_list_from_csv(args.input_file, args.csv_file)
    elif args.fingerprint_file:
        info = generate_list_from_fingerprint(args.input_file, args.fingerprint_file)

with open(args.output_file, 'a+') as report:
    report.write(info)
