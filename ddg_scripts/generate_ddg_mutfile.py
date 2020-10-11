import argparse
from Bio import SeqIO
import datetime
from os import mkdir, makedirs, remove, rename
from os.path import isdir, isfile, join
import pandas as pd
from pyrosetta import *
from pyrosetta.rosetta.core.pose import Pose, get_chain_from_chain_id, remove_nonprotein_residues, \
    get_chain_id_from_chain
from pyrosetta.rosetta.core.sequence import SWAligner, Sequence, SimpleScoringScheme
import re
import sys


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template_pdb', type=str, required=True, 
        help='Input a starting PDB file for comparison and from which mutants \
        will be generated.')
    parser.add_argument('-m', '--mutants_list', type=str, nargs='*', required=True, 
        help='Input a .fasta list file or files identifying the mutations. For an \
        oligomer protein, you may want to give multiple input .fasta files \
        corresponding to different chains in the pdb model. To do this, use space \
        to separate different .fasta files.')
    parser.add_argument('-cut', '--cut_region_by_chains', type=str, nargs='*', 
        help='if multiple fasta files input, cut regions are needed to be defined \
        in the same order of fasta files order. example: "A C B"')
    parser.add_argument('-dup', '--duplicated_chains', type=str, nargs='*', 
        help='Declare if the protein is a symmmetric protein.')
    return parser.parse_args()

def find(s, ch): # find all characters and return their indices
    return [i for i, ltr in enumerate(s) if ltr == ch]

def convert_date_str(date_str):
    """ 
    Converts a date in the form 'yyyy-mm-dd' to a datetime date object. Hyphens 
    can be any non-numeric separator.
    """
    # Convert date to integer list
    date_space = re.sub('[^0-9]', ' ', date_str) # Change separators to spaces
    date_list = [int(i) for i in date_space.split()]

    # Output date object
    return datetime.date(*date_list)

def read_name_tag(fasta_id):
    """
    Takes a biopython SeqRecord object with a fasta ID in the following form: 
    2020-03-28|2020-03-28|Count=1|hCoV-19/USA/WA-S424/2020|hCoV-19/USA/WA-S424/2020|NSP5
    Separates each |-separated feature to store as an attribute, with '/' 
    replaced by '_'. Returns a dict with all values.
    """
    tag_dict = {}

    breakup_id = fasta_id.split('|')
    print(breakup_id)
    tag_dict['date_first'] = convert_date_str(breakup_id[0])
    tag_dict['date_last'] = convert_date_str(breakup_id[1])
    tag_dict['count'] = int(breakup_id[2].split('=')[-1])
    tag_dict['id_1'] = breakup_id[3]
    tag_dict['id_2'] = breakup_id[4]
    tag_dict['location_1'] = breakup_id[3].split('/')[1]
    tag_dict['location_2'] = breakup_id[4].split('/')[1]
    tag_dict['tag'] = breakup_id[5]
    
    return tag_dict

def get_id(fasta_name):
    """
    Added by Changpeng. Used for check_replicates function
    """
    fp = open(fasta_name,"r")
    labels = []
    n = 0
    for line in fp:
        if line[0] == ">":
            if n != 0:
                labels.append(line.strip()[1:].split("|")[3]) # Find replicates by ID1
        n += 1
    fp.close()
    return labels

def parse_fastafile(fasta_file):
    """
    Read in a file with a list of fasta sequences and return a list of biopython
    SeqRecord objects for all sequences
    """
    # Initialize list
    fasta_list = []

    # Populate list
    for r in SeqIO.parse(fasta_file, 'fasta'): 
        fasta_list.append(r)

    return fasta_list

def replicate_seqs(replicates, analyze_lists):
    """
    Find the replicate sequences' indices in the fasta files. The output format is shown below:
    {{'hCoV-19/South_Africa/R05475/2020':SeqRecord1, SeqRecord2}
    """
    if replicates != None:
        id1 = []
        inds = {}
        for lt in analyze_lists:
            id1.append([str(x.id).split("|")[3] for x in lt])
        for rep in replicates.keys():
            inds[rep]=([analyze_lists[y][id1[y].index(rep)] for y in list(replicates[rep])])
    else:
        inds = None
    return inds

def seq_length_by_chain(wild_pose):
    """
    Added by Changpeng.
    Get the sequence length for each chain. Format is a dictionary, an example is shown below:
    ['DKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIV \
    QLSEISMDNSPNLAWPLIVTALRA',840], ['KMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQG', 954]]
    """
    lengths_chain_id = []
    wild_seqs = {} #"every value is: [seq, start_ind]"
    former_len = 0
    for chain in range(1,wild_pose.num_chains()+1):
        chain_name = get_chain_from_chain_id(chain, wild_pose)
        chain_seq = wild_pose.chain_sequence(chain)
        if chain_name in wild_seqs.keys():
            wild_seqs[chain_name][0] += wild_pose.chain_sequence(chain)
            #wild_seqs[chain_name][1] += len(chain_seq)
            #former_len += len(chain_seq)
        else:
            wild_seqs[chain_name] = []
            wild_seqs[chain_name].append(chain_seq)
            wild_seqs[chain_name].append(former_len + 1)
        former_len += len(chain_seq)
    return wild_seqs

def cut_by_chain(wild_pose, cut, list_fasta_names):
    """
    if multiple fasta files, then cut wild_pose into several parts that are matched with FASTA files.
    cut is based on cut_region_by_chains flag. Format is a dictionary, an example is shown below:
    ['DKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIV
    QLSEISMDNSPNLAWPLIVTALRA',840], ['KMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQG', 954]]
    """
    wild_seqs = {}
    tmp_pose = wild_pose.clone()
    remove_nonprotein_residues(tmp_pose)
    #if len(list_fasta_name) > 1:
    if cut:
        chain_seqs = seq_length_by_chain(wild_pose)
        try:
            if len(cut) == len(list_fasta_names):
                for re in cut:
                    chain = get_chain_id_from_chain(re, tmp_pose)
                    wild_seqs[re] = [tmp_pose.chain_sequence(chain), chain_seqs[re][1]]                   
        except TypeError:
            print("Invalid cut regions specified! Please add '-cut' flag or make sure the number \
            of regions are the same as the number of FASTA files! Be caution to place regions in order \
            with FASTA files after the '-m' flag")
    else:
        wild_seqs[get_chain_from_chain_id(1,wild_pose)] = [wild_pose.chain_sequence(1),1]
    return wild_seqs


def check_replicates(list_fasta_names):
    """
   Added by Changpeng. Find replicate reference proteins through multiple FASTA input files by searching for the same ID1.
   Output will be a dictionary that is the format like this:
   {id1: set(0,1), id1: set(0,1,2)}
   """
    if len(list_fasta_names) > 1:
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
                    preset.append((i,j))
                    preset.append((j,i))
                    indices.append((i,j))
        replicates = {}
        preset = []
        for mapping in indices:
            for que in ids[mapping[0]]:
                try:
                    tmp_ind = ids[mapping[1]].index(que)
                    if que not in preset:
                        replicates[que] = set([mapping[0],mapping[1]])
                        preset.append(que)
                    else:
                        replicates[que] = set(list(replicates[que]) + [mapping[0],mapping[1]])
                except ValueError:
                    pass
    else:
        replicates = None
    return replicates

def compare_sequences(pdb_name, pdb_seq, seq_2, query_seq, ind_by, is_False): #ind_by has two option, pose or pdb
    """
    Given a reference sequence and a comparison sequence, identify the sites 
    where the comparison sequence differs from the reference. Returns a list of 
    3-member tuples. Each trio represents one point substitution. The three  
    items listed in the trio are the site, the starting residue type, and the  
    substitute residue type. Site numbers are 1-indexed (not 0-indexed as Python 
    default) to work with Rosetta. If the sequences are not the same length, it  
    breaks the script.
    """
    # Make sure sequences are the same length
    # added by Changpeng, record the start index in pdb file
    if ind_by == "pdb":
        fp = open(pdb_name,"r")
        for line in fp:
            if line[0:6].strip() == "ATOM":
                if line[21] == list(pdb_seq.keys())[0]:
                    start_ind = int(line[22:26].strip()) # record the start index of residues
                    break
        fp.close()
    elif ind_by == "pose":
        start_ind = 1
    print("Start index in the pdb sequence is:{}".format(start_ind))

    # added by Changpeng, do sequence alignment between pdb seq and ref seq
    wild_seq = list(pdb_seq.values())[0][0] #one chain that matched to fasta,pdb_seq[1] is the start_ind in the wild_pose
    wild_seq_ = Sequence(wild_seq,"pdb")
    len_wild_seq = len(wild_seq)
    chain_name = list(pdb_seq.keys())[0]
    start_ind_chain_in_pose = list(pdb_seq.values())[0][1]
    
    seq_2 = str(seq_2).strip("-")
    seq_2_ = Sequence(seq_2, "reference")
    query_seq_ = Sequence(str(query_seq.seq).strip("-"), str(query_seq.id)) 
    ss = SimpleScoringScheme()
    tmpAlign = SWAligner().align(wild_seq_, seq_2_, ss)

    ss = SimpleScoringScheme()
    tmpAlign1 = SWAligner().align(wild_seq_, query_seq_, ss)

    align_ind_1 = int(tmpAlign.to_string().strip().split("\n ")[1].split()[1])
    align_ind_2 = int(tmpAlign.to_string().strip().split("\n ")[2].split()[1])
    new_seq_1 = tmpAlign.to_string().strip().split("\n ")[1].strip().split(" ")[-1].replace("\n", "")
    new_seq_2 = tmpAlign.to_string().strip().split("\n ")[2].strip().split(" ")[-1].replace("\n", "")
    query_seq_aligned = tmpAlign1.to_string().strip().split("\n")[2].strip().split(" ")[-1].replace("\n","")

    if is_False:
        align_ind_1 = 1
        align_ind_2 = start_ind_chain_in_pose
        new_seq_1 = ""
    #   new_seq_2 = seq_2[start_ind_chain_in_pose-1:]
        print(len(wild_seq))
        fp = open(pdb_name,"r")
        j = 0
        n = 0
        for line in fp:
            if line[0:6].strip() == "ATOM":
                if line[21] == list(pdb_seq.keys())[0]:
                    if int(line[22:26].strip()) > n:
                        if int(line[22:26].strip()) != n + 1:
                            new_seq_1 += "-" * (int(line[22:26].strip()) - n - 1)
                        new_seq_1 += wild_seq[j]
                        j += 1
                        n = int(line[22:26].strip())
        new_seq_1 = new_seq_1.lstrip("-")
        new_seq_2 = seq_2[start_ind-1: n] + "--"
        query_seq_aligned = str(query_seq.seq).strip("-")[start_ind-1: n] + "--"

    site_changes = []

    # N-terminal has truncations or additions in fasta
    if align_ind_1 != 1 or align_ind_2 != 1:
        if align_ind_1 < align_ind_2:
            for more_ind in range(0, align_ind_2 - align_ind_1):
                if query_seq[more_ind] != seq_2[more_ind]: 
                    site_changes.append((start_ind - (align_ind_2 - align_ind_1 ) + more_ind, query_seq[more_ind], seq_2[more_ind], False, chain_name)) # label that additions on N terminal for reference sequences; add one more column to judge if it is in pdb or not
        elif align_ind_1 > align_ind_2:
            new_seq_1 = wild_seq[0 : align_ind_1 - align_ind_2] + new_seq_1
            new_seq_2 = wild_seq[0 : align_ind_1 - align_ind_2] + new_seq_2
            query_seq_aligned = wild_seq[0 : align_ind_1 - align_ind_2] + query_seq_aligned

    # # C-terminal truncations has already been truncated in aligned seq_1, record C-terminal truncations
    if new_seq_2[-1] == "-":
        if len(new_seq_1.replace("-","")) != len(wild_seq): # need to be modified generally, len(new_seq_1) > len(new_seq_2)
            len_in_wild = len(new_seq_1.replace("-", ""))
            len_before_1 = len(new_seq_1)
            new_seq_1 = new_seq_1 + wild_seq[len_before_1:] # Now N terminal has been adjusted
            new_seq_2 = new_seq_2.strip("-") + wild_seq[len_before_1:]
            query_seq_aligned = query_seq_aligned.strip("-") + wild_seq[len_before_1:]

        else:
            new_seq_2 = new_seq_2.strip("-")
            query_seq_aligned = query_seq_aligned.strip("-")
            start_ind_in_seq2 = len(new_seq_2)
            if align_ind_1 < align_ind_2:
                start_ind_in_seq2 += align_ind_2 - 1
            if len(find(new_seq_2, "-")) != 0:
                start_ind_in_seq2 -= len(find(new_seq_2,"-"))
            for co_C in range(start_ind_in_seq2, len(seq_2)):
                if query_seq[co_C] != seq_2[co_C]:
                    site_changes.append((start_ind + co_C - (align_ind_2 - align_ind_1) + len(find(new_seq_2, "-")), query_seq[co_C], seq_2[co_C], False, chain_name)) # just to symbol len(ref) > len(wild)
    assert len(new_seq_1) == len(new_seq_2)
    print("PDB sequence after alignment is:")
    print(new_seq_1)
    print("Reference sequence after alignment is:")
    print(new_seq_2)
    print("Wild Type sequence after alignment is:")
    print(query_seq_aligned)
    # Record site changes that remove gaps for subsequent Rosetta repack + minimization
    new_site_changes = []
    # if gaps in reference sequences
    hyphen_list = find(new_seq_1, "-")
    seq_2_for_model = ''.join([new_seq_2[sel] for sel in range(len(new_seq_2)) if sel not in hyphen_list])        
    seq_1_for_model = new_seq_1.replace("-","")
    print("\n")
    print(seq_2_for_model)
    print(seq_1_for_model)
    hyphen_list2 = find(seq_2_for_model, "-")
    if len(hyphen_list2) != 0:
        for xx in hyphen_list2:
            seq_2_for_model = seq_2_for_model[0:xx] + seq_1_for_model[xx] + seq_2_for_model[xx+1:]
    print(len(seq_2_for_model))
    print(len(seq_1_for_model))

    assert len(seq_1_for_model) == len(seq_2_for_model)
    # Identify altered sites
    for n, i in enumerate(new_seq_2):
        if i != query_seq_aligned[n]:
            site = n + start_ind
            seq_1_value = query_seq_aligned[n]
            seq_2_value = i
            site_changes.append((site, seq_1_value, seq_2_value, new_seq_1[n] != "-", chain_name))

    for nn, ii in enumerate(seq_2_for_model):
        if ii != seq_1_for_model[nn]:
            site = nn + start_ind_chain_in_pose
            seq_1_value = seq_1_for_model[nn]
            seq_2_value = ii
            new_site_changes.append((site, seq_1_value, seq_2_value))
    print(site_changes)
    print(new_site_changes)
    return site_changes, new_site_changes

def generate_substitutuions(seqrecord, query, pdb_seq, fa_ind, pdb_name, cut_order, replicate_id1, rep_fa_ind, rep_searched):
    # Break up FASTA tag to component data
    mut_tags = read_name_tag(seqrecord.id)

    # Identify sequence differences from the original
    print("Chain is:{}".format(list(fa_ind.keys())[0]))
    print("This is the {}th FASTA file!".format(list(fa_ind.values())[0]))
    substitutions, new_subs = compare_sequences(pdb_name, {list(fa_ind.keys())[0] : pdb_seq[list(fa_ind.keys())[0]]}, seqrecord.seq, 
        query[list(fa_ind.values())[0]], 'pdb', False)

    # replicates processing
    if replicate_id1 != None:
        if mut_tags['id_1'] in replicate_id1.keys():
            fastas = replicate_id1[mut_tags['id_1']]
            fa_inds = list(rep_fa_ind[mut_tags['id_1']])
            fastas.pop(0)
            fa_inds.pop(0)
            for ff, rest in enumerate(fastas):
               print("Chain is:{}".format(cut_order[fa_inds[ff]]))
               print("This is the {}th FASTA file!".format(fa_inds[ff]))
               res_substitutions, res_new_subs = compare_sequences(pdb_name, {cut_order[fa_inds[ff]]: pdb_seq[cut_order[fa_inds[ff]]]}, \
                    rest.seq, query[fa_inds[ff]], 'pdb', False)
               substitutions += res_substitutions
               new_subs += res_new_subs
            replicate_id1.pop(mut_tags['id_1'])
    
    rep_searched.append(mut_tags['id_1'])
    print("For protein {0}, substitutions are {1}".format(mut_tags['id_1'],new_subs))
    return new_subs, rep_searched

if __name__ == '__main__':
    args = parse_arguments()

    init()
    wild_pose = pose_from_pdb(args.template_pdb)
    print("------------------------------")
    print("Considered FASTA files are:")
    print(args.mutants_list)

    replicates = check_replicates(args.mutants_list)

    #Read all FASTA files. need to combine with check_replicates together further
    analyze_lists = []
    wts = []
    for i in args.mutants_list:
        fasta_list = parse_fastafile(i)
        analyze_lists.append(fasta_list[1:])
        wts.append(fasta_list[0])
    # get the replicate inds in FASTA files
    replicate_inds = replicate_seqs(replicates,analyze_lists)
    print(replicate_inds)
    
    # Take subset of fasta list if parallelizing
    #analyze_list = partition_list(fasta_list[1:], *args.parallel_partition)
    #print(analyze_list[0])

    #wt = fasta_list[0]
    # split wild_pose into different parts that are corresponding to cut_regions
    pdb_seqs = cut_by_chain(wild_pose, args.cut_region_by_chains, args.mutants_list)
    print("Cut regions of the pdb and their start index in the pdb are:")
    print(pdb_seqs)

    if not args.cut_region_by_chains:
        args.cut_region_by_chains = list(pdb_seqs.keys())[0]

    # Iterate over identified fasta sequences, altering protease model
    all_mutants_info = pd.DataFrame([])
    all_substitutions = []
    replicates_searched = []
    protein_name = args.template_pdb[:args.template_pdb.find('_')]
    # iterate over fasta files
    for c in range(len(args.mutants_list)):
        fingerprint_file_list = list()
        total_mutations = 0
        mut_file_list = list()
        # iterate over sequences in current fasta file
        for n, mutant in enumerate(analyze_lists[c]):
            if read_name_tag(mutant.id)['id_1'] not in replicates_searched:
                pdb_seqs = cut_by_chain(wild_pose, args.cut_region_by_chains, args.mutants_list)
                new_subs, replicates_searched = generate_substitutuions(mutant, wts, \
                    pdb_seqs, {args.cut_region_by_chains[c]: c}, args.template_pdb, \
                    args.cut_region_by_chains, replicate_inds, replicates, replicates_searched)
                # Fingerprint for wild type or mutated variants
                fingerprint = str()
                # Only calculate mutated variants in Cartesian ddG
                if len(new_subs) > 0:
                    # Append all point mutations of the variant to mutfile
                    mut_list = list()
                    # iterate over point mutations
                    for pm in new_subs:
                        native_res = wild_pose.residue(pm[0]).name1()
                        if native_res == pm[1]:
                            fingerprint += pm[1] + str(pm[0]) + pm[2] + ';'
                            mut_list.append(pm[1] + ' ' + str(pm[0]) + ' ' + pm[2])
                        else:
                            print('The native residue type at pose position ' + \
                                str(pm[0]) + ' is ' + native_res + ' instead of ' + pm[1])
                            native_res_2 = wild_pose.residue(pm[0] + 1).name1()
                            if native_res_2 == pm[1]:
                                fingerprint += pm[1] + str(pm[0] + 1) + pm[2] + ';'
                                mut_list.append(pm[1] + ' ' + str(pm[0] + 1) + ' ' + pm[2])
                            else:
                                raise Exception('The native residue type at pose position ' + \
                                    str(pm[0]) + ' is ' + native_res + ' instead of ' + pm[1])
                        total_mutations += 1
                        # duplicate point mutations if has duplicated chains
                        if args.duplicated_chains:
                            res_pdb_info = list(filter(lambda x: x != '', wild_pose.pdb_info().\
                                pose2pdb(pm[0]).split(' ')))
                            # Current point mutation is not matched from other fasta files
                            if res_pdb_info[1] == args.duplicated_chains[0]:
                                for duplicated_chain in args.duplicated_chains[1:]:
                                    duplicated_point_mutation_pose_index = wild_pose.pdb_info().\
                                        pdb2pose(duplicated_chain, int(res_pdb_info[0]))
                                    fingerprint += pm[1] + str(duplicated_point_mutation_pose_index) + pm[2] + ';'
                                    mut_list.append(pm[1] + ' ' + str(duplicated_point_mutation_pose_index) + ' ' + pm[2])
                                    total_mutations += 1
                    # end of for loop for point mutations
                    mut_file_list.append(mut_list)
                else:
                    fingerprint += 'WT,'
                fingerprint_file_list.append(fingerprint[:-1] + '\n')
            else:
                print("Replicate that is already detected!. Skip this round of mutation.")
        # end of for loop for sequences in current fasta file
        # Single .fasta.txt file or multiple XXX_matched_0.fasta.txt files
        if len(args.mutants_list) > 1:
            last_idx = args.mutants_list[c].rfind('_')
            prefix = protein_name + '_' + args.mutants_list[c][last_idx + 1:-10]
        else: # Single XXX_matched_n.fasta.txt file, n > 0
            prefix = args.mutants_list[c][:-10]
        # fingerprint file would not be generated if all variants in current fasta file are matched to other fasta files
        if len(fingerprint_file_list) > 0:
            # Create the work directory
            mkdir(prefix)
            rename(args.mutants_list[c], prefix + '/' + args.mutants_list[c])
            # Write point mutation information of all variants to the fingerprint file
            p_fingerprint = open(prefix + '/' + prefix + '.fingerprint', 'a+')
            p_fingerprint.writelines(fingerprint_file_list)
            p_fingerprint.close()
        # mutfile would not be generated if all variants in current fasta file are wild type 
        # or all variants in current fasta file are matched to other fasta files
        if len(mut_file_list) > 0:
            # Write out to mutfile (here the mutant_list file name ends with the job index)
            if isfile(prefix + '/' + prefix + '.mut'):
                # Read the original mutfile
                p_mut = open(prefix + '/' + prefix + '.mut', 'r')
                mut_list_bak = p_mut.readlines()
                p_mut.close()
                # Remove the original mutfile
                remove(prefix + '/' + prefix + '.mut')
                # Calculate original total number of the point mutations
                total_mutations_old = int(mut_list_bak[0].split(' ')[1][:-1])
                # Rewrite the mutfile
                p_mut = open(prefix + '/' + prefix + '.mut', 'w+')
                # Calculate current total number of the point mutations
                total_mutations += total_mutations_old
                p_mut.write('total ' + str(total_mutations) + '\n')
                # Write the original part
                p_mut.writelines(mut_list_bak[1:])
                # Write current parts
                for mut_list in mut_file_list:
                    p_mut.write(str(len(mut_list)) + '\n')
                    for mut in mut_list:
                        p_mut.write(mut + '\n')
                p_mut.close()
            else:
                p_mut = open(prefix + '/' + prefix + '.mut', 'w+')
                p_mut.write('total ' + str(total_mutations) + '\n')
                for mut_list in mut_file_list:
                    p_mut.write(str(len(mut_list)) + '\n')
                    for mut in mut_list:
                        p_mut.write(mut + '\n')
                p_mut.close()
    # end of for loop for fasta files
