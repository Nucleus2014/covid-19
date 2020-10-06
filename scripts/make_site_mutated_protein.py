"""
Given a starting PDB structure and a list of mutant FASTA sequences perform 
several functions. 
Makes a table all substitution sites identified across the set:
    - The initial residue type
    - The identified substituent residues at that site
    - Whether the residue is in the core, boundary layer, or surface of the 
    protein
    - If there is a substrate in the model (assumed to be the second chain), 
    whether the substitution site is near it. 
    - If a catalytic residue is specified, whether the substitution site is 
    near it. 
Makes a table of all mutants: 
    - Input information about the collection of the mutant sequence 
    - Which substitution(s) is/are present
    - The energetic impact of the substitution(s)
    - Whether there are multiple substitutions
    - If there are multiple substitutions, whether they interact
Makes structural models of each substituted mutant
Optionally, the step of making the altered models can be omitted. If this option 
is chosen, the process will be considerably faster, skipping the most 
time-consuming steps. However, energy differences will not be calculated and the
potential interaction of mutations will be less accurately calculated from the 
wild-type model.
Script by Joseph H. Lubin, summer 2020
jhl133@scarletmail.rutgers.edu
"""

import argparse
from Bio import SeqIO
import datetime
import matplotlib.pyplot as plt
import numpy as np
from os import makedirs
from os.path import isdir, join
import pandas as pd
import pyrosetta as pr
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.pose import get_chain_from_chain_id, \
    get_chain_id_from_chain, Pose, remove_nonprotein_residues
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
    RestrictAbsentCanonicalAASRLT, RestrictToRepackingRLT, PreventRepackingRLT, \
    RestrictToRepacking
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector, \
    ChainSelector, LayerSelector, InterGroupInterfaceByVectorSelector, \
    NotResidueSelector, OrResidueSelector, ResidueIndexSelector, \
    TrueResidueSelector, NeighborhoodResidueSelector, ResiduePropertySelector
from pyrosetta.rosetta.core.sequence import SWAligner, Sequence, read_fasta_file,\
    SimpleScoringScheme
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    InteractionEnergyMetric, RMSDMetric, TotalEnergyMetric
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import \
    PerResidueEnergyMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.membrane import AddMembraneMover
from pyrosetta.rosetta.protocols.membrane.symmetry import \
    SymmetricAddMembraneMover
from pyrosetta.rosetta.protocols.minimization_packing import \
    MinMover, PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
import re
import sys

########## Collecting User Input ###############################################

def parse_args():
    info= """
    Given a starting PDB structure and a list of mutant FASTA sequences perform 
    several functions. 
    Makes a table all substitution sites identified across the set:
        - The initial residue type
        - The identified substituent residues at that site
        - Whether the residue is in the core, boundary layer, or surface of the 
        protein
        - If there is a substrate in the model (assumed to be the second chain), 
        whether the substitution site is near it. 
        - If a catalytic residue is specified, whether the substitution site is 
        near it. 
    Makes a table of all mutants: 
        - Input information about the collection of the mutant sequence 
        - Which substitution(s) is/are present
        - The energetic impact of the substitution(s)
        - Whether there are multiple substitutions
        - If there are multiple substitutions, whether they interact
    Makes structural models of each substituted mutant
    Optionally, the step of making the altered models can be omitted. If this 
    option is chosen, the process will be considerably faster, skipping the most 
    time-consuming steps. However, energy differences will not be calculated and 
    the potential interaction of mutations will be less accurately calculated by 
    distance instead of energy.
    """
    parser = argparse.ArgumentParser(description=info)
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
    parser.add_argument('-od', '--out_dir', type=str, 
        help='Input a directory into which the homolog models will be saved. \
        If not specified, PDBs will be saved in the current directory.')
    parser.add_argument('-rn', '--report_name', type=str, 
        help='Input a name for the substitutions report. \
        If not specified, the report will be called substitutions_summary.csv \
        and be saved in the current directory.')
    parser.add_argument('-params', '--params', type=str, nargs='+', 
        default=None, help='If a non-canonical residue/ligand is present, \
        provide a params file.')
    parser.add_argument('-mc', '--main_chain', type=int, default=1, 
        help='The main chain being analyzed is 1 by default. Specify main \
        chain if not 1.')
    parser.add_argument('-ic', '--interface_chain', type=int, nargs='+', 
        default=None, help='If symmetry is used, specify the oligomeric chains.')
    parser.add_argument('-lig', '--ligand_chain', type=int, nargs='+', 
        default=None, help='If a ligand is present, specify ligand chains. If \
        using both symmetry and ligands, specify only ligand chains that are \
        bound to the main chain.')
    parser.add_argument('-cr', '--catalytic_residues', type=int, nargs='+', 
        default=None, help='The catalytic residues of the enzyme. By default, \
        no residues are so designated. If residues are specified, report will \
        include whether substitutions interact with the catalytic residues.')
    parser.add_argument('-pmm', '--is_pdb_index_match_mutant', action='store_true', \
        help = 'if pdb index matches mutant index, then set this flag to avoid pairwise \
        alignment which costs much time and not so accurate')
    parser.add_argument('-ind', '--ind_type', choices = ['pose', 'pdb'], 
        default = 'pdb', help='To show mutation residues indices in pdb or in pose order')
    parser.add_argument('-sym', '--symmetry', type=str,
        help='If the pose is symmetric, include a symdef file.')
    parser.add_argument('-memb', '--membrane', required=False, action='store_true',
        help='Declare if the protein is a membrane protein.')
    parser.add_argument('-mspan', '--span_file', 
        required=any(x in ['--membrane','-memb'] for x in sys.argv),
        help='If the pose is a membrane protein, include a spanfile.')
    parser.add_argument('-proto', '--protocol', choices=['fastrelax', 'repack+min'], 
        default = 'repack+min', help='Employ either fast relax protocol or repacking \
        and minimization on the pose.')
    parser.add_argument('-ite', '--iterations', type=int, default=1, 
        help='Giving a flag of -ite, [ite] trajectories will be run for each variant \
        and output the decoy with the lowest score.')
    parser.add_argument('-rep', '--repulsive_type', type=str, nargs=2, choices=['soft', 'hard'], 
        default=['hard', 'hard'], help='Using the normal hard-rep or the \
        soft-rep score function for repacking and minimization, respectively.')
    parser.add_argument('-rnd', '--rounds', type=int, default=1, 
        help='The rounds of repacking and minimization calculations being repeated \
        in one trajectory.')
    parser.add_argument('-rep_wts', '--repulsive_weights', type=float, nargs='*', \
        default=None, help='Ramping Lennard-Jones repulsive weight in each rounds.')
    parser.add_argument('-no_ex', '--extra_rotamers', action='store_false', 
        help='Giving a flag of -no_ex will prevent using extra rotamer \
        sampling during repacking, saving a lot of time with reduced sampling.')
    parser.add_argument('-nbh', '--neighborhood_residue', type=float, 
        help='Giving a flag of -nbh will also allow surrounding residues \
            within [nbh] angstroms of the mutated residue to repack.')
    parser.add_argument('-op', '--only_protein', action='store_true', 
        help='Giving a flag of -op will prevent non-protein motifs from repacking, \
        such as ligands and RNA.')
    parser.add_argument('-fix_bb', '--backbone', action='store_false', 
        help='Giving a flag of -fix_bb will hold the backbone fixed in minimization.')
    parser.add_argument('-no_cst', '--constrain', action='store_false', 
        help='Giving a flag of -no_cst will prevent coordinate constraints \
        from being applied to the pose during repacking and minimization.')
    parser.add_argument('-coev', '--coveloution_method', type=str, 
        choices=['energy', 'distance'], default='distance', 
        help='Whether mutations are a coevolution can be detected in one of two \
        ways: geometrically, based on the inter-residue distance and \
        orientation, or energetically, based on whether residues have nonzero \
        interface scores. Uses geometry by default.')
    parser.add_argument('-no_cst_score', '--unconstrain_ddg', action='store_true', \
        default=False, help="Call this flag to return energies unconstrained")
    parser.add_argument('-no_models', '--make_models', action='store_false', 
        help='Giving a flag of -no_models will prevent PDB models from being \
        generated. This will prevent the energetic calculations and not yield \
        athe substituted models, but will skip the most time-consuming steps.')
    parser.add_argument('-part', '--parallel_partition', type=int, nargs=2, 
        default=[1,1], help='To parallelize the analysis, enter two integers. \
        The first is the number of partitions, the second is which partition \
        member to run on this processor, from 1 to the number of partitions. \
        not working for now.')
    parser.add_argument('-debug', '--debugging_mode', action='store_true', 
        help='Giving a flag of -debug, it will print out the task operations on \
        all residues, namely point mutations, repacking or keeping static.')
    return parser.parse_args()

########## Text-based data collection ##########################################
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
    print("Now processing: {}".format(fasta_id))
    breakup_id = fasta_id.split('|')
    tag_dict['date_first'] = convert_date_str(breakup_id[0])
    tag_dict['date_last'] = convert_date_str(breakup_id[1])
    tag_dict['gisaid_count'] = int(breakup_id[2].split('=')[-1])
    tag_dict['id_1'] = breakup_id[3]
    tag_dict['id_2'] = breakup_id[4]
    tag_dict['location_1'] = breakup_id[3].split('/')[1]
    tag_dict['location_2'] = breakup_id[4].split('/')[1]
    tag_dict['tag'] = breakup_id[5]
    
    return tag_dict


def compare_sequences(pdb_name, pdb_seq, seq_2, query_seq, ind_by, 
    is_pdb_index_match_mutant): 
    """
    Given a reference sequence and a comparison sequence, identify the sites 
    where the comparison sequence differs from the reference. Returns a list of 
    3-member tuples. Each trio represents one point substitution. The three  
    items listed in the trio are the site, the starting residue type, and the  
    substitute residue type. Site numbers are 1-indexed (not 0-indexed as Python 
    default) to work with Rosetta. If the sequences are not the same length, it  
    breaks the script.

    ind_by has two options, pose or pdb
    """
    # Make sure sequences are the same length
    # added by Changpeng, record the start index in pdb file
    if ind_by == "pdb":
        fp = open(pdb_name,"r")
        for line in fp:
            if line[0:6].strip() == "ATOM":
                if line[21] == list(pdb_seq.keys())[0]:
                    # record the start index of residues
                    start_ind = int(line[22:26].strip()) 
                    break
        fp.close()
    elif ind_by == "pose":
        start_ind = 1
    print("Start index in the pdb sequence is:{}".format(start_ind))

    # added by Changpeng, do sequence alignment between pdb seq and ref seq
    #one chain that matched to fasta,pdb_seq[1] is the start_ind in the wild_pose
    wild_seq = list(pdb_seq.values())[0][0] 
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

    if is_pdb_index_match_mutant:
        align_ind_1 = 1
        align_ind_2 = start_ind_chain_in_pose
        new_seq_1 = ""
    #    new_seq_2 = seq_2[start_ind_chain_in_pose-1:]
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


def get_mutant_brief(point_substitution):
    """
    Input a 5-member tuple representing a point substitution, of the form: 
    (site, original_AA, substituted_AA, whether_site_is_in_pdb, chain)
    Note: this is the output of the compare_sequences function.

    Returns the typical string summary of a point substitution, such as 
    A_S328A, which indicates that site 328 on chain A mutated from S to A.
    """
    template = '{}_{}{}{}'
    ps = point_substitution
    m = template.format(ps[4], ps[1], ps[0], ps[2])

    return m


def is_conservative(point_substitution):
    """
    Input a 3-member tuple representing a point substitution, of the form: 
    (site, original_AA, substituted_AA)
    Note: this is the output of the compare_sequences function.
    Returns a Boolean of whether a change from the first residue to the second
    represents a conservative substitution.
    Uses groupings of amino acids from Wikipedia:
    https://en.wikipedia.org/wiki/Amino_acid#/media/File:Amino_Acids.svg
    Staying within B, or D will be conservative. 
    Staying withing charge group in A will be conservative. 
    Any special case is nonconservative.
    """
    res_groups = [  ['A', 'F', 'I', 'L', 'M', 'V', 'W', 'Y'],
                    ['C'],
                    ['D', 'E'],
                    ['G'], 
                    ['H','K','R'],
                    ['N', 'Q', 'S', 'T'],
                    ['P'],
                    ['U']]

    # Identify in which group the original residue belongs
    for rg in res_groups:
        if point_substitution[1] in rg:
            original_group = rg
            break

    # Check whether the substituted residue is in the same group
    conservation = point_substitution[2] in original_group

    return conservation

########## Rosetta-based data collection #######################################

def selector_intersection(*selectors):
    """ Returns the intersection of any set of selectors """
    intersect_selection = AndResidueSelector()
    for s in selectors:
        intersect_selection.add_residue_selector(s)

    return intersect_selection


def selector_union(*selectors):
    """ Returns the intersection of any set of selectors """
    union_selection = OrResidueSelector()
    for s in selectors:
        union_selection.add_residue_selector(s)

    return union_selection


def selector_to_list(pose, selector):
    """
    Produces a list of residues from a pose that are identified 
    by a given selector
    """
    # Make selection
    selection = selector.apply(pose)
    
    # Convert to a vector of selected residues
    selected_vector = get_residues_from_subset(selection)
    
    # Convert vector to a list
    selected_list = list(selected_vector)
    
    return selected_list


def identify_res_layer(pose, res_number, main_chain=1):
    """
    Determines whether a given residue in a pose is in the core, boundary, or 
    surface layer of the protein.
    """
    # If the PDB has multiple chains, isolate main
    check_pose = pr.Pose(pose, pose.chain_begin(main_chain), pose.chain_end(main_chain))
    
    # Identify layer with LayerSelector
    layer_selector = LayerSelector()

    # Checking core
    layer_selector.set_layers(1, 0, 0)
    core_selection = layer_selector.apply(pose)
    if core_selection[res_number]:
        return 'core'

    # Checking boundary
    layer_selector.set_layers(0, 1, 0)
    boundary_selection = layer_selector.apply(pose)
    if boundary_selection[res_number]:
        return 'boundary'

    # Checking surface
    layer_selector.set_layers(0, 0, 1)
    surface_selection = layer_selector.apply(pose)
    if surface_selection[res_number]:
        return 'surface'


def total_energy(pose, score_function, selection=None):
    """
    Calculates total energy of a pose using a TotalEnergyMetric. If a selector
    is provided, calculates the total energy of the selection rather than the
    whole pose.
    """
    # Create the metric
    tem = TotalEnergyMetric()
    tem.set_scorefunction(score_function)
    
    # Add the selector
    if selection:
        tem.set_residue_selector(selection)
        
    return tem.calculate(pose)


def get_rmsd(pose_1, pose_2):
    """
    Given two poses of equal size, determines RMSD. The first pose is the 
    one to which the second is compared.
    """
    # Create the RMSDMetric, setting pose_1 as the reference
    rmsd_metric = RMSDMetric()
    rmsd_metric.set_comparison_pose(pose_1)

    # Use the RMSDMetirc to calculate the RMSD of pose_2
    rmsd = rmsd_metric.calculate(pose_2)
    
    return rmsd


def interaction_energy(pose, score_function, selection_1, selection_2):
    """
    Given a Rosetta pose, score function, and two residue selectors, 
    calculates the interaction energy between the selections.
    """
    interact_metric = InteractionEnergyMetric()
    interact_metric.set_scorefunction(score_function)
    interact_metric.set_residue_selectors(selection_1, selection_2)
    return interact_metric.calculate(pose)


def check_selection_proximity(pose, selection_1, selection_2):
    """
    Determines whether any residues in one selector are close enough to 
    potentially interact with residues from a second selector, 
    using an InterGroupInterfaceByVector selector with default settings.
    """
    # Make full-pose selection
    full_pose = TrueResidueSelector()
    
    # Make selection for the full pose minus the first selector
    not_selection_1 = NotResidueSelector(selection_1)
    full_minus_1 = selector_intersection(full_pose, not_selection_1)
    
    # Make selection for interacting residues between selection 1 the rest of the pose
    igibv = InterGroupInterfaceByVectorSelector()
    igibv.group1_selector(selection_1)
    igibv.group2_selector(full_minus_1)
    
    # Find intersection of selector_1's interaction partners and selector_2
    selection_overlap = selector_intersection(igibv, selection_2)
    
    # Check if there are residues in the intersection, or if it is empty
    # If it is empty, the Boolean will be False. Otherwise, it will be True
    selections_do_overlap = bool(selector_to_list(pose, selection_overlap))

    return selections_do_overlap


def check_coevolution(pose, substitutions, score_function=None):
    """
    Given a pose and and a list of substitution sites, makes selectors for the 
    substitution sites and checkes whether the selections interact. If a score
    function is given, interaction is defined by having a nonzero interface 
    energy. If no score function is given, interaction is defined as having CA
    within 5.5A or within 11A if the CBs are oriented toward each other.
    """
    # Set output to False unless an interaction is found.
    sites_interact = False

    # Automatically no interaction if there are not multiple substitutions
    if len(substitutions) < 2:
        return sites_interact

    # Non-redundantly iterate pairwise through the list of substitutions 
    # Skip last res to avoid self-comparison
    for n, sub1 in enumerate(substitutions[:-2]):
        # Make selector for first residue 
        res1 = ResidueIndexSelector(str(sub1))

        # Iterate through all subsequent substitution sites
        for sub2 in substitutions[n + 1:]:
            # Make selector for second residue 
            res2 = ResidueIndexSelector(str(sub2))

            # Checking interaction by energy
            if score_function:
                # Determine interface energy
                intE = interaction_energy(pose, score_function, res1, res2)

                # If the energy is nonzero, residues are interacting
                if intE != 0:
                    sites_interact = True
                    break 

            # Checking interaction by geometry
            else:
                # Calculating proximity
                selections_nearby = check_selection_proximity(pose, res1, res2)

                # If selections are nearby, residues are interacting
                if selections_nearby:
                    sites_interact = True
                    break 


    # When either an interaction is found or when iteration completes, return
    # whether sites interact
    return sites_interact


def check_interface_proximity(pose, score_function, 
    substitute_res, sym_chains, main_chain=1):
    """
    Given a pose and a score function, and an integer for substitute_res
    and a list of chains for sym_chains, makes selectors for the substituted 
    residue and for the symmetric chains, and also for the main chain. Generates 
    an interface selection between the main chain and symmetric chains, then 
    checkes whether the substituted residue and interface selections have 
    nonzero interaction energy, indicating that the substitution is in the
    interfacial region between symmetric chains.
    """
    # Make a selector for the substituted residue
    mutated_res_selection = ResidueIndexSelector(substitute_res)

    # Make a selector for the symmetric chains
    individual_chain_selections = [ChainSelector(i) for i in sym_chains]
    sym_chain_selections = selector_union(*individual_chain_selections)

    # Make selector for the main chain
    main_chain_selection = ChainSelector(main_chain)

    # Making a selection for the interface residues between chains
    oligomer_interface = InterGroupInterfaceByVectorSelector()
    oligomer_interface.group1_selector(main_chain_selection)
    oligomer_interface.group2_selector(sym_chain_selections)

    # Limiting interface selection only to residues off the main chain
    not_main = NotResidueSelector(main_chain_selection)
    partner_interface = selector_intersection(oligomer_interface, not_main)

    # Check interaction energy between selections
    interchain_interaction_energy = interaction_energy(pose, score_function, 
        mutated_res_selection, partner_interface)
    interface_proximity_by_energy = bool(interchain_interaction_energy)

    return interface_proximity_by_energy


def check_catalytic_proximity(pose, score_function, 
    substitute_res, catalytic_res):
    """
    Given a pose and a score function, and an integer for substitute_res
    and a list of integers for catalytic_res, makes selectors for the 
    substitute_res and for the catalytic_res, and checkes whether the selections
    have nonzero interaction energy.
    """
    # Make a selector for the substituted residue
    mutated_res_selection = ResidueIndexSelector(substitute_res)

    # Make a selector for the catalytic residues
    cat_string = ','.join([str(i) for i in catalytic_res])
    catalytic_selection = ResidueIndexSelector(cat_string)

    # Check interaction energy between selections
    catalytic_interaction_energy = interaction_energy(pose, score_function, 
        mutated_res_selection, catalytic_selection)
    catalytic_proximity_by_energy = bool(catalytic_interaction_energy)

    return catalytic_proximity_by_energy

########## Making Substituted Models ########################################### 

def coord_constrain_pose(pose):
    """
    Applies full coordinate constraints to a pose
    Modified by Zhuofan according to Elliott's suggestion
    """
    cg = CoordinateConstraintGenerator()
    cg.set_bounded_width(0.1)
    cg.set_bounded(True)
    cg.set_sidechain(False)
    cg.set_sd(0.5)
    ac = AddConstraints()
    ac.add_generator(cg)
    ac.apply(pose)

    return pose


def make_point_mutant_task_factory(input_pose, site_changes, ex12=True, repacking_range=False, 
    only_protein=False, duplicated_chains=None):
    """
    Input a site_changes list of 3-member tuples, each representing a point 
    substitution in the form output by the compare_sequences function: 
    (site, original_AA, substituted_AA)
    Generates a TaskFactory to repack all residues in a pose with the 
    altered sites repacking to the new sequence.
    The ex12 option perform extra sampling ar chi 1 and chi 2 when 
    repacking. This increased sampling substantially increases runtime.
    Modified by Zhuofan
    """
    # Make TaskFactory to input changes
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    if ex12:
        tf.push_back(ExtraRotamers(0, 1, 1))
        tf.push_back(ExtraRotamers(0, 2, 1))

    # Force packing to new residue for each mutated site
    mutated_res_selection = OrResidueSelector() # Initialize selection of all altered residues
    for pm in site_changes:
        res_selection = ResidueIndexSelector(str(pm[0]))
        aa_force = RestrictAbsentCanonicalAASRLT()
        aa_force.aas_to_keep(pm[2])
        tf.push_back(OperateOnResidueSubset(aa_force, res_selection))
        mutated_res_selection.add_residue_selector(res_selection)
        # duplicate point mutations if has duplicated chains
        if duplicated_chains:
            res_pdb_info = list(filter(lambda x: x != '', input_pose.pdb_info().\
                pose2pdb(pm[0]).split(' ')))
            # Current point mutation is not matched from other fasta files
            if res_pdb_info[1] == duplicated_chains[0]:
                for duplicated_chain in duplicated_chains[1:]:
                    duplicated_point_mutation_pose_index = input_pose.pdb_info().\
                        pdb2pose(duplicated_chain, int(res_pdb_info[0]))
                    res_selection = ResidueIndexSelector(duplicated_point_mutation_pose_index)
                    aa_force = RestrictAbsentCanonicalAASRLT()
                    aa_force.aas_to_keep(pm[2])
                    tf.push_back(OperateOnResidueSubset(aa_force, res_selection))
                    mutated_res_selection.add_residue_selector(res_selection)
    # Repack some residues to accommodate substitutions
    repack = RestrictToRepackingRLT()
    prevent = PreventRepackingRLT()
    if repacking_range:
        # Repack surrounding residues within [repacking_range] angstroms
        repacking_res_selection = NeighborhoodResidueSelector()
        repacking_res_selection.set_focus_selector(mutated_res_selection)
        repacking_res_selection.set_distance(repacking_range)
        repacking_res_selection.set_include_focus_in_subset(False)
        if only_protein:
            protein_selection = ResiduePropertySelector(ResidueProperty.PROTEIN)
            repacking_res_selection_only_protein = AndResidueSelector(repacking_res_selection, \
                protein_selection)
            tf.push_back(OperateOnResidueSubset(repack, repacking_res_selection_only_protein))
            mutated_and_repacking_res_selection = OrResidueSelector(mutated_res_selection, \
                repacking_res_selection_only_protein)
            tf.push_back(OperateOnResidueSubset(prevent, mutated_and_repacking_res_selection, True))
        else:
            tf.push_back(OperateOnResidueSubset(repack, repacking_res_selection))
            # Do nothing to other residues
            mutated_and_repacking_res_selection = OrResidueSelector(mutated_res_selection, repacking_res_selection)
            # By setting the arg2 to True, the mutation and repacking res selection is flipped
            tf.push_back(OperateOnResidueSubset(prevent, mutated_and_repacking_res_selection, True))
    else:
        if only_protein:
            protein_selection = ResiduePropertySelector(ResidueProperty.PROTEIN)
            repack_res_selection = AndResidueSelector(NotResidueSelector(mutated_res_selection), protein_selection)
            tf.push_back(OperateOnResidueSubset(repack, repack_res_selection))
            tf.push_back(OperateOnResidueSubset(prevent, protein_selection, True))
        else:
            # By setting the arg2 to True, the mutation selection is flipped
            tf.push_back(OperateOnResidueSubset(repack, mutated_res_selection, True))

    return tf


def make_ref_pose_task_factory(input_pose, site_changes, ex12=True, repacking_range=False, 
    only_protein=False, duplicated_chains=None):
    """
    Written by Zhuofan
    """
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    if ex12:
        tf.push_back(ExtraRotamers(0, 1, 1))
        tf.push_back(ExtraRotamers(0, 2, 1))
    repack = RestrictToRepackingRLT()
    prevent = PreventRepackingRLT()
    if repacking_range:
        mutated_res_selection = ResidueIndexSelector()
        for pm in site_changes:
            mutated_res_selection.append_index(int(pm[0]))
            # duplicate point mutations if has duplicated chains
            if duplicated_chains:
                res_pdb_info = list(filter(lambda x: x != '', input_pose.pdb_info().\
                    pose2pdb(pm[0]).split(' ')))
                # Current point mutation is not matched from other fasta files
                if res_pdb_info[1] == duplicated_chains[0]:
                    for duplicated_chain in duplicated_chains[1:]:
                        duplicated_point_mutation_pose_index = input_pose.pdb_info().\
                            pdb2pose(duplicated_chain, int(res_pdb_info[0]))
                        mutated_res_selection.append_index(duplicated_point_mutation_pose_index)
        repacking_res_selection = NeighborhoodResidueSelector()
        repacking_res_selection.set_focus_selector(mutated_res_selection)
        repacking_res_selection.set_distance(repacking_range)
        repacking_res_selection.set_include_focus_in_subset(True)
        if only_protein:
            protein_selection = ResiduePropertySelector(ResidueProperty.PROTEIN)
            repacking_res_selection_only_protein = AndResidueSelector(repacking_res_selection, \
                protein_selection)
            tf.push_back(OperateOnResidueSubset(repack, repacking_res_selection_only_protein))
            tf.push_back(OperateOnResidueSubset(prevent, repacking_res_selection_only_protein, True))
        else:
            tf.push_back(OperateOnResidueSubset(repack, repacking_res_selection))
            tf.push_back(OperateOnResidueSubset(prevent, repacking_res_selection, True))
    else:
        if only_protein:
            protein_selection = ResiduePropertySelector(ResidueProperty.PROTEIN)
            tf.push_back(OperateOnResidueSubset(repack, protein_selection))
            tf.push_back(OperateOnResidueSubset(prevent, protein_selection, True))
        else:
            tf.push_back(RestrictToRepacking())

    return tf


def make_move_map(pose, backbone, only_protein):
    """
    Generates a movemap for either (Minimization after Repacking) or 
    (only FastRelax protein scaffold).
    Written by Zhuofan
    """
    move_map = pr.MoveMap()
    move_map.set_bb(backbone)
    if only_protein:
        protein_selection = ResiduePropertySelector(ResidueProperty.PROTEIN)
        protein_res_true_vector = protein_selection.apply(pose)
        move_map.set_chi(protein_res_true_vector)
    else:
        move_map.set_chi(True)
    move_map.set_jump(True)
    return move_map


def repacking_with_muts_and_minimization(pose, score_functions, decoys, rounds, 
    repulsive_weights, is_ref_pose, site_changes, ex12=True, repacking_range=False, 
    only_protein=False, duplicated_chains=None, backbone=True):
    """
    Applies point mutations to a given pose. This is done through a 
    PackRotamersMover, followed by minimization.
    pose is a Rosetta pose
    site_changes is a list of point_mutation objects
    Modified by Zhuofan.
    """
    # Set up a repacker
    pack_rotamer_mover = PackRotamersMover()
    # Set up a flexible backbone or fixed backbone minimizer
    min_mover = MinMover()
    min_mover.min_type('lbfgs_armijo_nonmonotone')
    move_map = make_move_map(pose, backbone, only_protein)
    min_mover.movemap(move_map)
    min_mover.min_type('lbfgs_armijo_nonmonotone')
    
    # Repeating the whole protocol for n decoys and get the best decoy.
    for decoy in range(decoys):
        # Make a copy Pose
        pose_copy = pr.Pose(pose)
        
        # Create a task factory for each copy of the pose object
        if is_ref_pose: # Make a task factory for reference pose
            task_factory = make_ref_pose_task_factory(pose_copy, site_changes, \
                ex12=ex12, repacking_range=repacking_range, only_protein=only_protein, \
                duplicated_chains=duplicated_chains)
        else: # Make a task factory to substitute residues and repack
            task_factory = make_point_mutant_task_factory(pose_copy, site_changes, \
                ex12=ex12, repacking_range=repacking_range, only_protein=only_protein, \
                duplicated_chains=duplicated_chains)
        # print(task_factory.create_task_and_apply_taskoperations(pose))
        # Set the task factory for the repacker.
        pack_rotamer_mover.task_factory(task_factory)

        # Repeating repacking and minimization for n rounds.
        for rnd in range(rounds):
            # Set score functions for the repacker and minimizer.
            if repulsive_weights: # Set fa_rep weight.
                repacking_score_function = score_functions[0].clone()
                repacking_score_function.set_weight(ScoreType.fa_rep, repulsive_weights[rnd])
                pack_rotamer_mover.score_function(repacking_score_function)
                minimization_score_function = score_functions[1].clone()
                minimization_score_function.set_weight(ScoreType.fa_rep, repulsive_weights[rnd])
                min_mover.score_function(minimization_score_function)
            else: # Use the default fa_rep weight.
                pack_rotamer_mover.score_function(score_functions[0])
                min_mover.score_function(score_functions[1])
            
            # Apply the PackRotamersMover
            pack_rotamer_mover.apply(pose_copy)
            # Apply the MinMover
            min_mover.apply(pose_copy)
        
        # Save the best decoy.
        current_energy = total_energy(pose_copy, score_functions[1])
        if decoy == 0 or current_energy < lowest_energy:
            lowest_energy = current_energy
            mutated_pose = pose_copy

    return mutated_pose


def fast_relax_with_muts(pose, score_function, decoys, is_ref_pose, site_changes, 
    ex12=True, repacking_range=False, only_protein=False, duplicated_chains=None, 
    backbone=True):
    """
    Applies point mutations to a given pose. This is done through a 
    Fast Relax Mover.
    Written by Zhuofan
    """
    # Set up a FastRelax mover. Set repeating trajectories to 1.
    fast_relax = FastRelax(1)
    fast_relax.set_scorefxn(score_function)
    move_map = make_move_map(pose, backbone, only_protein)
    fast_relax.set_movemap(move_map)

    # Repeating the FastRelax protocol for n decoys and get the best decoy.
    for decoy in range(decoys):
        # Make a copy Pose
        pose_copy = pr.Pose(pose)
        
        # Create a task factory for each copy of the pose object
        if is_ref_pose: # Make a task factory for reference pose
            task_factory = make_ref_pose_task_factory(pose_copy, site_changes, \
                ex12=ex12, repacking_range=repacking_range, only_protein=only_protein, \
                duplicated_chains=duplicated_chains)
        else: # Make a task factory to substitute residues and repack
            task_factory = make_point_mutant_task_factory(pose_copy, site_changes, \
                ex12=ex12, repacking_range=repacking_range, only_protein=only_protein, \
                duplicated_chains=duplicated_chains)
        # print(task_factory.create_task_and_apply_taskoperations(pose))
        # Set the task factory for the FastRelax mover.
        fast_relax.set_task_factory(task_factory)

        # Apply the FastRelax mover.
        fast_relax.apply(pose_copy)

        # Save the best decoy.
        current_energy = total_energy(pose_copy, score_function)
        if decoy == 0 or current_energy < lowest_energy:
            lowest_energy = current_energy
            mutated_pose = pose_copy

    return mutated_pose


def make_mutant_model(ref_pose, substitutions, score_functions, 
    protocol='repack+min', decoys=1, rounds=1, repulsive_weights=None, 
    ex12=True, repacking_range=False, backbone=True, only_protein=False, 
    no_constraint_scoring=True, duplicated_chains=None):
    """
    Given a reference pose and a list of substitutions, and a score function, 
    produces a mutated pose that has the substitutions and is repacked and 
    minimized. Also outputs a dict with the energy change and RMSD of the 
    mutated pose
    substitutions should be list of 3-member tuples, each representing a point 
    substitution in the form output by the compare_sequences function: 
    (site, original_AA, substituted_AA)
    The ex12 option will perform extra sampling ar chi 1 and chi 2 when 
    repacking. This increased sampling substantially increases runtime.
    Modified by Zhuofan.
    """
    if len(substitutions) == 0: # No change.
        mutated_pose = pr.Pose(ref_pose)
        modified_ref_pose = pr.Pose(ref_pose)
    else: # Make residue changes.
        # Double check the point mutation informations
        corrections_to_substitutions = dict()
        for idx, pm in enumerate(substitutions):
            native_res = ref_pose.residue(pm[0]).name1()
            if native_res != pm[1]:
                raise Exception('The native residue type at pose position ' + \
                    str(pm[0]) + ' is ' + native_res + ' instead of ' + pm[1])
        # Possibly make corrections to new_sub
        for idx, new_point_mutation_info_tuple in corrections_to_substitutions.items():
            substitutions[int(idx)] = new_point_mutation_info_tuple
        
        if protocol == 'repack+min':
            mutated_pose = repacking_with_muts_and_minimization(ref_pose, \
                score_functions, decoys, rounds, repulsive_weights, True, \
                substitutions, ex12=ex12, repacking_range=repacking_range, \
                only_protein=only_protein, duplicated_chains=duplicated_chains, \
                backbone=backbone)
            modified_ref_pose = repacking_with_muts_and_minimization(ref_pose, \
                score_functions, 1, rounds, repulsive_weights, False, \
                substitutions, ex12=ex12, repacking_range=repacking_range, \
                only_protein=only_protein, duplicated_chains=duplicated_chains, \
                backbone=backbone)
        elif protocol == 'fastrelax':
            mutated_pose = fast_relax_with_muts(ref_pose, score_functions[1], decoys, \
                True, substitutions, ex12=ex12, repacking_range=repacking_range, \
                only_protein=only_protein, duplicated_chains=duplicated_chains, \
                backbone=backbone)
            modified_ref_pose = fast_relax_with_muts(ref_pose, score_functions[1], 1, \
                False, substitutions, ex12=ex12, repacking_range=repacking_range, \
                only_protein=only_protein, duplicated_chains=duplicated_chains, \
                backbone=backbone)

    # Switch back to a single score function
    score_function = score_functions[1]

    # Initialize data collection dict
    mutated_pose_data = {}

    # Rescore the poses
    wt_energy = score_function(modified_ref_pose)
    mutant_energy = score_function(mutated_pose)

    # Remove constraint weights from total energy
    if no_constraint_scoring:
        cst_type = ScoreType.coordinate_constraint
        coord_wt = score_function.get_weight(ScoreType.coordinate_constraint)
        wild_cst = modified_ref_pose.energies().total_energies().get(cst_type)
        wt_energy = wt_energy - wild_cst * coord_wt
        mutant_cst = mutated_pose.energies().total_energies().get(cst_type)
        mutant_energy = mutant_energy - mutant_cst * coord_wt
    
    # Add energy change to output data
    mutated_pose_data['energy_change'] = mutant_energy - wt_energy

    # Add RMSD to output data
    mutated_pose_data['rmsd'] = get_rmsd(modified_ref_pose, mutated_pose)

    return modified_ref_pose, mutated_pose, mutated_pose_data

########## Collecting Model Data ###############################################

def out_directory(directory_name):
    """
    Given an output directory string, checks whether that directory exists.  
    Creates it if it doesn't. If no name is given, returns an empty string.
    """
    if directory_name:
        outdir = directory_name
        if not isdir(directory_name):
            makedirs(directory_name)
    else:
        outdir = ''

    return outdir


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


def partition_list(in_list, partitions, member):
    """
    Given a list, a number of partitions, and a partition member number,
    splits the list into the appropriate number of partitions, and returns 
    the specified member. Member should be 1-indexed, so the minimum is 1 and
    maximum is the number of partitions.
    """
    # Confirm appropriate member number
    assert 1 <= member <= partitions

    # Determine list length and partition size
    list_size = len(in_list)
    partition_size = int(list_size/partitions)
    overrun = list_size % partitions

    # Determine starting index to collect. If the list doesn't break into equal 
    # partitions, the earlier partitions will have one extra element
    start_index = 0
    for i in range(1, member):
        if i <= overrun:
            start_index += partition_size + 1
        else:
            start_index += partition_size

    # Determine end index to collect
    if member <= overrun:
        end_index = start_index + partition_size + 1
    else:
        end_index = start_index + partition_size

    return in_list[start_index:end_index]


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
        chain_seq = Pose(wild_pose, wild_pose.chain_begin(chain), wild_pose.chain_end(chain)).sequence()
        if chain_name in wild_seqs.keys():
            wild_seqs[chain_name][0] += Pose(wild_pose, wild_pose.chain_begin(chain), wild_pose.chain_end(chain)).sequence()
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
        print(chain_seqs)
        try:
            if len(cut) == len(list_fasta_names):
                for re in cut:
                    chain = get_chain_id_from_chain(re, tmp_pose)
                    align_ind = 0
                    true_start_ind = chain_seqs[re][1]
                    pseudo_chain_seq = chain_seqs[re][0]
                    while pseudo_chain_seq[align_ind] != tmp_pose.residue(tmp_pose.chain_begin(chain)).name1():
                        true_start_ind += 1
                        pseudo_chain_seq  = pseudo_chain_seq[1:]
                    print(pseudo_chain_seq)
                    wild_seqs[re] = [Pose(tmp_pose, tmp_pose.chain_begin(chain), \
                        tmp_pose.chain_end(chain)).sequence(), true_start_ind]
        except TypeError:
            print("Invalid cut regions specified! Please add '-cut' flag or make sure the number \
            of regions are the same as the number of FASTA files! Be caution to place regions in order \
            with FASTA files after the '-m' flag")
    else:
        wild_seqs[get_chain_from_chain_id(1,wild_pose)] = [wild_pose.chain_sequence(1),1]
    return wild_seqs


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


def analyze_mutant_protein(seqrecord, ref_pose, score_functions, query, pdb_seq, 
    fa_ind, pdb_name, make_model=True, main_chain=1, oligo_chains=None, 
    substrate_chains=None, cat_res=None, ind_type='pdb', 
    pdb_index_match_mutant=False, cut_order=None, rep_fa_ind=None, 
    replicate_id1=None, rep_searched=None, protocol='repack+min', decoys=1, 
    rounds=1, ex12=True, repacking_range=False, backbone=True, 
    only_protein=False, debugging_mode=False, unconstrain_final_score=True, 
    duplicated_chains=None, coveloution_method='distance'):
    """
    Given a biopython SeqRecord object with a fasta ID in the following form: 
    2020-03-28|2020-03-28|Count=1|hCoV-19/USA/WA-S424/2020|hCoV-19/USA/WA-S424/2020|NSP5
    The ID is then broken up so each |-separated feature is entered separately 
    into a dict. Also reports countries where the sequences were first and most 
    recently collected, what substitutions the sequence contains, whether 
    multiple substitutions are present, and whether substitutions are 
    conservative. 
    If reference pose and score function are provided, will report several 
    other checks. If make_model is true, a mutated model of the reference pose
    which includes the identified substitutions will be created, and all 
    evaluations will be performed on the mutated model.
    Other checks include:
    --If a mutated model is made, the energy change and RMSD compared to the 
        original model
    --The burial layer (Core, Boundary, Surface) of substituted residues
    --If multiple substitutions are present, whether they energetically interact
    --If symmetry chains are specified as a list in oligo_chains, checks whether 
        substitutions interact energetically with interfacial residues
    --If catalytic residues are specified as a list in cat_res, checks whether 
        substitutions interact energetically with catalytic residues
    --If ligands are specified as a list in substrate_chains, checks whether 
        substitutions interact energetically with ligands
    The ex12 indicates whether to perform extra sampling ar chi 1 and chi 2 when 
    repacking. This increased sampling substantially increases runtime. 
    """
    # Break up FASTA tag to component data
    mut_tags = read_name_tag(seqrecord.id)

    # Identify sequence differences from the original
    print("Chain is:{}".format(list(fa_ind.keys())[0]))
    print("This is the {}th FASTA file!".format(list(fa_ind.values())[0]))
    substitutions, new_subs = compare_sequences(pdb_name, 
        {list(fa_ind.keys())[0] : pdb_seq[list(fa_ind.keys())[0]]}, seqrecord.seq, 
        query[list(fa_ind.values())[0]], ind_type, pdb_index_match_mutant)

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
               res_substitutions, res_new_subs = compare_sequences(pdb_name, 
                    {cut_order[fa_inds[ff]]: pdb_seq[cut_order[fa_inds[ff]]]}, 
                    rest.seq, query[fa_inds[ff]], ind_type,pdb_index_match_mutant)
               substitutions += res_substitutions
               new_subs += res_new_subs
            replicate_id1.pop(mut_tags['id_1'])
 
    rep_searched.append(mut_tags['id_1'])
    print("For protein {0}, substitutions are {1}".format(mut_tags['id_1'],new_subs))
    
    # Add substitutions to output data
    sub_shorts = [get_mutant_brief(pm) for pm in substitutions]
    mut_tags['substitutions'] = ';'.join(sub_shorts)

    # Add whether mutation is in pose or not
    mut_tags['is_in_PDB'] = ';'.join([str(pm[3]).upper() for pm in substitutions])

    # Add to output data whether there are multiple substitutions
    mut_tags['multiple'] = len(new_subs) > 1

    # Add whether substitution(s) are conservative to output data
    conservations = [is_conservative(pm) for pm in substitutions]
    mut_tags['conservative'] = ';'.join([str(i).upper() for i in conservations])

    # If the mutant model is being made, the subsequent steps should use the 
    # Mutant model for increased accuracy, so that step is performed here
    if make_model:
        modified_ref_pose, mutated_pose, mut_pose_info = \
            make_mutant_model(ref_pose, new_subs, score_functions, 
                protocol=protocol, decoys=decoys, rounds=rounds, ex12=ex12, 
                repacking_range=repacking_range, backbone=backbone, 
                only_protein=only_protein, debugging_mode=debugging_mode,
                no_constraint_scoring=unconstrain_final_score, 
                duplicated_chains=duplicated_chains)
        
        # Switch back to a single score function
        sf = score_functions[1]

        # Add mutant model info to the output data
        mut_tags.update(mut_pose_info)

        # Check if is_in_pdb, replace NA with energy terms
        tmp = mut_tags['is_in_PDB'].split(";")
        if 'True' not in tmp:
            mut_tags['energy_change'] = 'NA'
            mut_tags['rmsd'] = 'NA'

        # Set future calculations to use the mutated pose
        ref_pose_copy = pr.Pose(ref_pose)
    else:
        mutated_pose = None
        ref_pose_copy = pr.Pose(ref_pose)
    # Add substitution layer to output data
    layers = [identify_res_layer(ref_pose_copy, pm[0], main_chain=list(fa_ind.values())[0]+1) 
        for pm in new_subs]
    mut_tags['layer'] = ';'.join([i.upper() for i in layers])

    # Add whether multiple substitutions are interacting to output data
    sites = [s[0] for s in new_subs]
    if coveloution_method == 'distance':
        coevolution = check_coevolution(ref_pose_copy, sites)

    if coveloution_method == 'energy':
        coevolution = check_coevolution(ref_pose_copy, sites, score_function=sf)
    mut_tags['coevolution'] = coevolution

    # Add whether substitutions are near multi-chain interface to output data
    if oligo_chains:
        interface_checks = [check_interface_proximity(ref_pose_copy, sf, pm[0], 
            sym_chains, main_chain=main_chain) for pm in new_subs]
        mut_tags['near_interface'] = ';'.join([str(i) for i in interface_checks])

    # Add whether substitutions are near catalytic residues to output data
    if cat_res:
        cat_checks = [check_catalytic_proximity(ref_pose_copy, sf, pm[0], cat_res)
            for pm in new_subs]
        mut_tags['near_catalytic'] = ';'.join([str(i) for i in cat_checks])

    # Add whether substitutions are near substrate to output data
    if substrate_chains:
        substrate_checks = [check_interface_proximity(ref_pose_copy, sf, pm[0], 
            substrate_chains, main_chain=main_chain) for pm in new_subs]
        mut_tags['near_substrate'] = ';'.join([str(i) for i in substrate_checks])

    # Convert to a DataFrame
    mut_df = pd.DataFrame(mut_tags, index=[1])

    return modified_ref_pose, mut_df, mutated_pose, substitutions, new_subs, replicate_id1, rep_searched

########## Main ################################################################         

def get_sf(rep_type, symmetry=False, membrane=False, constrain=True, cst_wt=1.0):
    """
    Determines the appropriate score function to use, based on a rep_type
    that is either hard (ref2015) or soft (ref2015_soft), whether symmetry 
    and/or membrane modeling are in use, and whether constraints are desired.
    """
    assert rep_type in ['hard', 'soft']

    if symmetry: # Declare symmetric score functions
        score_function = SymmetricScoreFunction()
        if rep_type == 'hard':
            if membrane:
                score_function.add_weights_from_file('franklin2019')
            else:
                score_function.add_weights_from_file('ref2015')
        elif rep_type == 'soft':
            if membrane: # Set up a soft-rep version of franklin2019 manually
                score_function.add_weights_from_file('ref2015_soft')
                score_function.set_weight(ScoreType.fa_water_to_bilayer, 1.0)
            else:
                score_function.add_weights_from_file('ref2015_soft')
    else: # Declare ordinary score functions
        if rep_type == 'hard':
            if membrane:
                score_function = pr.create_score_function('franklin2019')
            else:
                score_function = pr.create_score_function('ref2015')
        elif rep_type == 'soft':
            if membrane: # Set up a soft-rep version of franklin2019 manually
                score_function = pr.create_score_function('ref2015_soft')
                score_function.set_weight(ScoreType.fa_water_to_bilayer, 1.0)
            else:
                score_function = pr.create_score_function('ref2015_soft')
    # The score functions do not have constraint weights incorporated in 
    # themselves. If requisite, the constraint weights are added.
    if constrain:
        score_function.set_weight(ScoreType.atom_pair_constraint, cst_wt)
        score_function.set_weight(ScoreType.coordinate_constraint, cst_wt)
        score_function.set_weight(ScoreType.angle_constraint, cst_wt)
        score_function.set_weight(ScoreType.dihedral_constraint, cst_wt)
        score_function.set_weight(ScoreType.metalbinding_constraint, cst_wt)
        score_function.set_weight(ScoreType.chainbreak, cst_wt)
        score_function.set_weight(ScoreType.res_type_constraint, cst_wt)
        
    return score_function


def get_pose(pdb, symmetry=None, membrane=None, constrain=True):
    """
    Returns a pose from a PDB. 
    If a symmetry file is given, applies symmetry. 
    If a span file is given, applies membrane.
    """
    pose = pr.pose_from_file(pdb)
    
    # Applying symmetry if specified
    if symmetry: 
        sfsm = SetupForSymmetryMover(symmetry)
        sfsm.apply(pose)
            
    # Set up membrane for membrane protein
    if membrane:
        if symmetry:
            add_memb = pr.rosetta.protocols.membrane.symmetry.SymmetricAddMembraneMover(membrane)
            add_memb.apply(pose)
        else:
            add_memb = pr.rosetta.protocols.membrane.AddMembraneMover(membrane)
            add_memb.apply(pose)
            
    # Apply constraints
    if constrain:
        pose = coord_constrain_pose(pose)
    
    return pose


def print_pose_scores(pose):
    """
    Prints all energies in a pose that have a non-zero weight in REF2015
    Note that scores are unweighted.
    """
    nonzeros = ['fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_rep', 
        'fa_intra_sol_xover4', 'lk_ball_wtd', 'fa_elec', 'pro_close', 
        'hbond_sr_bb', 'hbond_lr_bb', 'hbond_bb_sc', 'hbond_sc', 'dslf_fa13', 
        'atom_pair_constraint', 'coordinate_constraint', 'angle_constraint', 
        'dihedral_constraint', 'metalbinding_constraint', 'omega', 'fa_dun', 
        'p_aa_pp', 'yhh_planarity', 'ref', 'chainbreak', 'rama_prepro', 
        'res_type_constraint', 'total_score']

    pose_energies = pose.energies().total_energies()
    listed_scores = pose_energies.show_nonzero().split()
    listed_scores = [i.rstrip(':') for i in listed_scores]
    for i in nonzeros:
        if i in listed_scores:
            print("{:<40}{}".format(i, listed_scores[listed_scores.index(i)+1]))
        else:
            print("{:<40}{}".format(i, 0))


def calc_single_res_Es(pose, sf, single_chain=True):
    """
    For a given pose and score function, returns a list of single-residue
    energies. If single_chain is true, only returns energies for the first 
    chain in the pose.
    """
    prem = PerResidueEnergyMetric()
    prem.set_scorefunction(sf)
    
    e_vals = prem.calculate(pose)

    if single_chain:
        return [e_vals[i] for i in range(1, pose.chain_end(1) + 1)]   
    else:
        return [e_vals[i] for i in range(1, pose.total_residue() + 1)]

def plot_res_ddG(ref_pose, mut_pose, sf, name=None):
    """
    Generates two per-residue energy plots from a pair of poses. The first
    is the total energies of each residue in both poses. The second is the 
    difference in per-residue energies between the poses. Reference is in red
    in the first plot, and the mut_pose is in blue.
    """
    ref_Es = calc_single_res_Es(ref_pose, sf)
    mut_Es = calc_single_res_Es(mut_pose, sf)
    dE = np.array(mut_Es) - np.array(ref_Es)

    # Plot
    fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(7.5, 5)
    fig.set_dpi(300)
    ax1.plot(ref_Es, linewidth=0.3, color='red')
    ax1.plot(mut_Es, linewidth=0.3, color='blue')
    ax2.plot(dE, linewidth=0.3, color='black')
    fig.tight_layout()
    if name:
        plt.savefig(name)
    else:
        plt.show()
    plt.close()

########## Executing ###########################################################

def main(args):
    # Load Rosetta
    opts = '-mute all -run:preserve_header'
    if args.params:
        opts += ' -extra_res_fa ' + ' '.join(args.params)
    # So far we don't allow switching the score function in fast relax
    if args.protocol == 'fastrelax':
        if args.repulsive_type:
            args.repulsive_type = None
        if args.repulsive_weights:
            args.repulsive_weights = None
    if args.protocol == 'repack+min':
        opts += ' -fa_max_dis 9.0'
    pr.init(opts)

    # Set up output directory if generating PDB models
    outdir = out_directory(args.out_dir)

    score_functions = list()
    for rep_type in args.repulsive_type:
        score_function = get_sf(rep_type, symmetry=args.symmetry, 
            membrane=args.membrane, constrain=args.constrain)
        score_functions.append(score_function)

    # Load protein model
    wild_pose = get_pose(args.template_pdb, symmetry=args.symmetry, 
        membrane=args.span_file, constrain=args.constrain)

    # Read fasta sequences
    # edited by Changpeng
    # Find replicates of reference proteins in multiple FASTA files
    # And read wild_seq for each monomer
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
        
    # split wild_pose into different parts that are corresponding to cut_regions
    pdb_seqs = cut_by_chain(wild_pose, args.cut_region_by_chains, args.mutants_list)
    print("Cut regions of the pdb and their start index in the pdb are:")
    print(pdb_seqs)

    if not args.cut_region_by_chains:
        args.cut_region_by_chains = list(pdb_seqs.keys())[0]

    # Iterate through all identified fasta sequences, altering protease model
    all_mutants_info = pd.DataFrame([])
    all_substitutions = []
    replicates_searched = []
    for c in range(len(args.mutants_list)):
        for n, mutant in enumerate(analyze_lists[c]):
            if read_name_tag(mutant.id)['id_1'] not in replicates_searched: 
                ref_pose = pr.Pose(wild_pose)
                modified_ref_pose, single_mutant_info, mutated_pose, \
                    substitutions, new_subs, replicate_inds, \
                    replicates_searched = analyze_mutant_protein(
                        mutant, ref_pose, score_functions, 
                        wts, pdb_seqs, {args.cut_region_by_chains[c]: c}, 
                        pdb_name=args.template_pdb, 
                        make_model=args.make_models, 
                        main_chain=args.main_chain, 
                        oligo_chains=args.interface_chain, 
                        substrate_chains=args.ligand_chain, 
                        cat_res=args.catalytic_residues, 
                        ind_type=args.ind_type, 
                        pdb_index_match_mutant=args.is_pdb_index_match_mutant, 
                        cut_order=args.cut_region_by_chains, 
                        rep_fa_ind=replicates, 
                        replicate_id1=replicate_inds, 
                        rep_searched=replicates_searched, 
                        protocol=args.protocol, 
                        decoys=args.iterations, 
                        rounds=args.rounds, 
                        ex12=args.extra_rotamers, 
                        repacking_range=args.neighborhood_residue, 
                        backbone=args.backbone, 
                        only_protein=args.only_protein, 
                        debugging_mode=args.debugging_mode, 
                        unconstrain_final_score=args.unconstrain_ddg, 
                        duplicated_chains=args.duplicated_chains, 
                        coveloution_method=args.coveloution_method)
            else:
                print("Replicate that is already detected!. Skip this round of mutation.")
            # Display results for this mutant
            if n == 0 and c == 0:
                print(single_mutant_info.to_string(header=True))
            else:
                print(single_mutant_info.to_string(header=False))

            # Output modified PDB
            if args.make_models:
                tag = single_mutant_info['tag'].values[0]
                is_pdb = single_mutant_info['is_in_PDB']
                first_id = single_mutant_info['id_1'].values[0].replace('/', '_')
                subs_clean = []
                for sub in substitutions:
                    #if sub[3] == True:
                        subs_clean.append((sub[4],sub[0],sub[1],sub[2]))
                mut_summary = [get_mutant_brief(s) for s in substitutions]
                str_mut_sum = '-'.join(mut_summary)
                outname = '{}_{}_{}.pdb'.format(tag, first_id, str_mut_sum)
                outname = join(outdir, outname)
                modified_ref_pose.dump_pdb(outname[:-4] + '_ref.pdb')
                mutated_pose.dump_pdb(outname)

                # Generate a per-residue energy comparison plot if debugging
                if args.debugging_mode:
                    plot_res_ddG(modified_ref_pose, mutated_pose,  
                        score_functions[1], outname[:-4] + \
                        '_res_score_changes.png')

            # Append individual mutant data to aggregate
            all_mutants_info = all_mutants_info.append(single_mutant_info)
            for s in subs_clean:
    #            if s not in all_substitutions:
                    all_substitutions.append(s)

    # Make report names
    if args.report_name:
        main_report_name = join(outdir, args.report_name + '_mutants')
        subs_report_name = join(outdir, args.report_name + '_substitutions')

    else:
        main_report_name = join(outdir, 'mutants_summary')
        subs_report_name = join(outdir, 'substitutions_summary')

    if args.parallel_partition != [1, 1]:
        main_report_name += '_{1}_of_{0}.csv'.format(*(args.parallel_partition))
        subs_report_name += '_{1}_of_{0}.csv'.format(*(args.parallel_partition))
    else:
        main_report_name += '.csv'
        subs_report_name += '.csv'

    # Make mutants summary csv
    all_mutants_info.to_csv(main_report_name, index=False)

    # Make substitutions summary a DataFrame, sort it, and output to csv
    counts = []
    all_subs_clean = list(set(all_substitutions))
    for sub in all_subs_clean:
        counts.append(np.sum(np.asarray([sub == x for x in all_substitutions],dtype=bool)))
    all_subs_info = pd.DataFrame(all_subs_clean)
    all_subs_info.columns = ['chain', 'site', 'native', 'mutant']
    all_subs_info['count'] = counts
    all_subs_info = all_subs_info.sort_values(by=['chain', 'site', 'mutant'])
    all_subs_info.to_csv(subs_report_name, index=False)


if __name__ == '__main__':
    args = parse_args()
    main(args)
