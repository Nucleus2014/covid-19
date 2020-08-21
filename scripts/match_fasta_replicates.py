import os
import argparse


def read_fasta(fasta_file_name, variant_dict, variant_list_list):
    #variant_dict = {"variant_name_3": [[fasta_name_a, line1, ...], [fasta_name_b line1, ...]], 
    #    "variant_name_2": [[fasta_name_a, line1, ...]], "variant_name_4": [[fasta_name_a, line1, ...]], 
    #    "variant_name_1": [[fasta_name_b, line1, ...]], ...}

    #seq_line_list_list = [[fasta_name_a, line1, ...], [fasta_name_b line1, ...]]

    #seq_line_list = [fasta_name_a, line1, ...]

    #variant_list = [variant_name_1, variant_name_2, variant_name_3, ...]

    with open(fasta_file_name, 'r') as p_fasta:
        # prefix of the fasta file name
        fasta_name = fasta_file_name[:-10]

        # Recording the original order of the variants in the fasta file
        variant_list = list()

        flag = 0
        for line in p_fasta:
            if line[0] == '>': # At the first line of a block

                if flag == 0: # At the first line of the reference bolck

                    flag = 1
                
                elif flag == 1: # At the first line of the first variant block

                    # Write the reference block to the matched_0 and matched fasta files.
                    with open(fasta_name +  '_matched_0.fasta.txt', 'w+') as new_fasta:
                        new_fasta.writelines(seq_line_list[1:])
                    with open(fasta_name + '_matched.fasta.txt', 'w+') as new_fasta:
                        new_fasta.writelines(seq_line_list[1:])
                    
                    # Get the name of this variant from the first line of this block
                    variant_name = line[1:].split('|')[3]
                    
                    variant_list.append(variant_name)

                    flag = 2
                
                elif flag == 2: # At the first line of other variant blocks

                    # Append the sequence lines list of the previous block to variant_dict[variant_name], 
                    # which is a list of the sequence lines list.
                    if not variant_dict.get(variant_name):
                        variant_dict[variant_name] = list()
                    seq_line_list_list = variant_dict[variant_name]
                    seq_line_list_list.append(seq_line_list)

                    # Get the name of this variant from the first line of this block
                    variant_name = line[1:].split('|')[3]

                    variant_list.append(variant_name)
                
                # if flag == 0 or flag == 1 or flag == 2:
                # After reading the first line of each block, create a sequence lines list
                seq_line_list = [fasta_name, line]
            
            else: # NOt at the first line of a block
                seq_line_list.append(line) # Append sequence lines to current sequence lines list 
        
        # For the very last block, there is no subsequent line starts with '>', so manually 
        # append the sequence lines list of the last block to variant_dict[variant_name].
        if not variant_dict.get(variant_name):
            variant_dict[variant_name] = list()
        seq_line_list_list = variant_dict[variant_name]
        seq_line_list_list.append(seq_line_list)

        # Save the list of variant names in this fasta file into 'fasta_variants_dict' 
        variant_list_list.append(variant_list)

def write_sequences(variant_dict, variant_list_list):
    #variant_dict = {"variant_name_3": [[fasta_name_a, line1, ...], [fasta_name_b line1, ...]], 
    #    "variant_name_2": [[fasta_name_a, line1, ...]], "variant_name_4": [[fasta_name_a, line1, ...]], 
    #    "variant_name_1": [[fasta_name_b, line1, ...]], ...}

    #seq_line_list_list = [[fasta_name_a, line1, ...], [fasta_name_b line1, ...]]

    #seq_line_list = [fasta_name_a, line1, ...]

    #variant_list = [variant_name_1, variant_name_2, variant_name_3, ...]

    for variant_list in variant_list_list:
        for variant in variant_list:
            seq_line_list_list = variant_dict.pop(variant, None)
            if seq_line_list_list:
                if len(seq_line_list_list) == 1:
                    seq_line_list = seq_line_list_list[0]
                    with open(seq_line_list[0] + '_matched.fasta.txt', 'a+') as p_fasta:
                        p_fasta.writelines(seq_line_list[1:])
                else:
                    for seq_line_list in seq_line_list_list:
                        with open(seq_line_list[0] + '_matched_0.fasta.txt', 'a+') as p_fasta:
                            p_fasta.writelines(seq_line_list[1:])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', type=str, nargs='*', help="List all fasta files that needs to be matched.")
    args = parser.parse_args()
    
    variants = dict()
    variants_mapping = list()
    for fasta_file in args.inputs:
        read_fasta(fasta_file, variants, variants_mapping)
    write_sequences(variants, variants_mapping)
