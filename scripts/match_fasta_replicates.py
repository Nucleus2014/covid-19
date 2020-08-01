import os
import argparse


def read_fasta(fasta_name, variant_dict):
    with open(fasta_name, 'r') as p_fasta:
        fasta_name = fasta_name[:-10]
        flag = 0
        for line in p_fasta:
            if line[0] == '>':
                if flag == 2:
                    if not variant_dict.get(variant_name):
                        variant_dict[variant_name] = list()
                    variant_dict[variant_name].append(seq_lines)
                    variant_name = line[1:].split('|')[3]
                elif flag == 1:
                    with open(fasta_name +  '_matched_0.fasta.txt', 'w+') as new_fasta:
                        new_fasta.writelines(seq_lines[1:])
                    with open(fasta_name + '_matched.fasta.txt', 'w+') as new_fasta:
                        new_fasta.writelines(seq_lines[1:])
                    variant_name = line[1:].split('|')[3]
                    flag = 2
                else:
                    flag = 1
                seq_lines = [fasta_name, line]
            else:
                seq_lines.append(line)
    if not variant_dict.get(variant_name):
        variant_dict[variant_name] = list()
    variant_dict[variant_name].append(seq_lines)

def write_sequences(variant_dict):
    for variant in variant_dict.values():
        for motif in variant:
            if len(variant) == 1:
                fasta_name = motif[0] + '_matched.fasta.txt'
            else:
                fasta_name = motif[0] + '_matched_0.fasta.txt'
            with open(fasta_name, 'a+') as p_fasta:
                p_fasta.writelines(motif[1:])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', type=str, nargs='*', help="List all fasta files that needs to be matched.")
    args = parser.parse_args()
    variants = dict()
    for fasta in args.inputs:
        read_fasta(fasta, variants)
    write_sequences(variants)
