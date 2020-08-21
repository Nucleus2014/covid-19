import argparse
from Bio import SeqUtils
from pyrosetta import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template_pdb', type=str, required=True, 
        help='Input a starting PDB file for comparison and from which mutants \
        will be generated.')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, 
        help='The output .csv file from make_site_mutated_protein.py.')
    parser.add_argument('-ddgs', '--ddg_files', type=str, nargs='*',
        help='The output .ddg files from previous calculations.')
    parser.add_argument('-n', '--mutfile_length', type=int, default=10, 
        help='Number of variants assigned to each .mut file.')
    parser.add_argument('-o', '--output_prefix', type=str, 
        help='Prefix of output .mut files.')
    args = parser.parse_args()
    if not args.output_prefix:
        args.output_prefix = args.template_pdb[:-4]
    return args

def read_variants_from_csv(csv):
    with open(csv, 'r') as p2csv:
        lines = p2csv.readlines()
    for idx, item in enumerate(lines[0].split(',')):
        if item == 'substitutions':
            break
    variants = list() # ['AI559_V;AP585_S', 'AT434_I'] (using pdb index)
    for line in lines[1:]:
        variant = line.split(',')[idx]
        if variant != '':
            variants.append(variant)
    return variants # ['AI559_V;AP585_S', 'AT434_I'] (using pdb index)

def read_variants_from_ddg(ddgs):
    ddg_variant_set = set()
    for ddg in ddgs:
        with open(ddg, 'r') as p2ddg:
            lines = p2ddg.readlines()
        for i in range(3, len(lines), 3):
            info = list(filter(lambda x: x != '', lines[i].split(' ')))
            ddg_variant_set.add(info[2][:-1])
    return ddg_variant_set

def convert_pdb_to_pose_numbering(pose, variants, ddg_variant_set=None): # variants always do not contain empty strings
    mut_list_list = list() # [['A 559 V\n', 'P 585 S\n'], ['T 434 I\n']] (using pose index)
    for variant in variants:
        mut_list = list() # ['A 559 V\n', 'P 585 S\n'] (using pose index)
        mut_str = 'MUT' # 'MUT_559VAL_585SER'
        mutations = variant.split(';')
        # For each point mutation, if it is inside the pdb, convert pdb index to pose index and append to mut_list and mut_str
        for mutation in mutations:
            chain_id = mutation[0]
            native_aa = mutation[1]
            mutated_aa = mutation[-1]
            pdb_index = int(mutation[2:-2])
            pose_index = pose.pdb_info().pdb2pose(chain_id, pdb_index)
            if pose_index != 0:
                mut_list.append(native_aa + ' ' + str(pose_index) + ' ' + mutated_aa + '\n')
                mut_str += '_' + str(pose_index) + SeqUtils.IUPACData.protein_letters_1to3[mutated_aa].upper()
        if mut_str == 'MUT': # pseudo wild type variant (point mutations locates in the pdb)
            print('============' + variant + ' is a pseudo wild type variant============\n')
            continue
        if ddg_variant_set:
            if mut_str in ddg_variant_set: # calculated in previous calculation (optional)
                print('============' + mut_str + ' (pose index) has already been calculated============\n')
                continue
        mut_list_list.append(mut_list)
    if ddg_variant_set:
        print('============Total number of non-wild type variants (not including pseudo wild-type variants) has been calculated is ' + str(len(ddg_variant_set)) + '============\n')
    print('============Total number of non-wild type variants (including pseudo wild-type variants) being analyzed is ' + str(len(variants)) + '============\n')
    print('============Total number of variants still need to be calculated is ' + str(len(mut_list_list)) + '============\n')
    return mut_list_list # [['A 559 V\n', 'P 585 S\n'], ['T 434 I\n']] (using pose index)

def write_mutfiles(mut_list_list, output_prefix, mutfile_length=1):
    # mut_list_list = [['A 559 V\n', 'P 585 S\n'], ['T 434 I\n']] (using pose index)
    mutfile_num = int(len(mut_list_list) / mutfile_length)
    for i in range(mutfile_num):
        mutfile_name = output_prefix + '_' + str(i + 1) + '.mut'
        begin = i * mutfile_length
        end = (i + 1) * mutfile_length
        write_mutfile(mutfile_name, mut_list_list, begin, end)
    last_mutfile_length = len(mut_list_list) % mutfile_length
    if last_mutfile_length > 0:
        mutfile_name = output_prefix + '_' + str(mutfile_num + 1) + '.mut'
        begin = mutfile_num * mutfile_length
        end = len(mut_list_list)

def write_mutfile(mutfile_name, mut_list_list, begin, end):
    # count total mutations of all variants
    total_mutations = 0
    for mut_list in mut_list_list[begin:end]:
        total_mutations += len(mut_list)
    # write each mutfile
    with open(mutfile_name, 'w+') as p2mut:
        # write total mutations of all variants
        p2mut.write('total ' + str(total_mutations) + '\n')
        for mut_list in mut_list_list[begin:end]:
            # count total mutations of each variant
            p2mut.write(str(len(mut_list)) + '\n')
            # write down point mutations
            p2mut.writelines(mut_list)    


if __name__ == '__main__':
    args = parse_args()
    init()
    pose = pose_from_pdb(args.template_pdb)
    variants = read_variants_from_csv(args.csv_file)
    if args.ddg_files:
        ddg_variant_set = read_variants_from_ddg(args.ddg_files)
        mut_list_list = convert_pdb_to_pose_numbering(pose, variants, ddg_variant_set)
    else:
        mut_list_list = convert_pdb_to_pose_numbering(pose, variants)
    write_mutfiles(mut_list_list, args.output_prefix, args.mutfile_length)
