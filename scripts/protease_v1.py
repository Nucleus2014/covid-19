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
import sys
from Bio import SeqIO
import datetime 
from os import makedirs
from os.path import isdir, join
import pandas as pd
import pyrosetta as pr
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.pose import get_chain_from_chain_id, remove_nonprotein_residues, \
    get_chain_id_from_chain
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
    RestrictAbsentCanonicalAASRLT, RestrictToRepackingRLT, PreventRepackingRLT
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
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.minimization_packing import \
    MinMover, PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
import re

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
    parser.add_argument('-m', '--mutants_list', type=str, required=True, 
        help='Input a fasta list file or files identifying the mutations.')
    parser.add_argument('-nbh', '--neighborhood_residue', type=float, 
        help='Giving a flag of -nbh will also allow surrounding residues \
            within [nbh] angstroms of the mutated residue to repack.')
    parser.add_argument('-fr', '--fast_relax', type=int, 
        help='Giving a flag of -fr will employ fast relax protocol on whole \
            protein instead of repacking and minimization, running for [fr] \
            trajectories.')
    parser.add_argument('-od', '--out_dir', type=str, 
        help='Input a directory into which the homolog models will be saved. \
        If not specified, PDBs will be saved in the current directory.')
    parser.add_argument('-rn', '--report_name', type=str, 
        help='Input a name for the substitutions report. \
        If not specified, the report will be called substitutions_summary.csv \
        and be saved in the current directory.')
    parser.add_argument('-mc', '--main_chain', type=int, default=1, 
        help='The main chain being analyzed is 1 by default. Specify main \
        chain if not 1.')
    parser.add_argument('-ic', '--interface_chain', type=int, nargs='+', 
        default=None, help='If symmetry is used, specify the oligomeric chains')
    parser.add_argument('-lig', '--ligand_chain', type=int, nargs='+', 
        default=None, help='If a ligand is present, specify ligand chains. If \
        using both symmetry and ligands, specify only ligand chains that are \
        bound to the main chain.')
    parser.add_argument('-sym', '--symmetry', type=str,
        help='If the pose is symmetric, include a symdef file.')
    parser.add_argument('-memb', '--membrane', required=False, action='store_true',
        help='Declare if the protein is a membrane protein.')
    parser.add_argument('-mspan', '--span_file', required='--membrane' or '-memb' in sys.argv,
        help='If the pose is a membrane protein, include a spanfile.')
    parser.add_argument('-cr', '--catalytic_residues', type=int, nargs='+', 
        default=None, help='The catalytic residues of the enzyme. By default, \
        no residues are so designated. If residues are specified, report will \
        include whether substitutions interact with the catalytic residues.')
    parser.add_argument('-params', '--params', type=str, nargs='+', 
        default=None, help='If a non-canonical residue/ligand is present, \
        provide a params file.')
    parser.add_argument('-cut', '--cut_region_by_chains', type=str, default=None, \
        help='if multiple fasta files input, cut regions are needed to be defined \
        in the same order of fasta files order. example: "A,C,B"')
    parser.add_argument('-no_cst', '--constrain', action='store_false', 
        help='Giving a flag of -no_cst will prevent coordinate constraints \
        from being applied to the pose during repacking and minimization.')
    parser.add_argument('-no_ex', '--extra_rotamers', action='store_false', 
        help='Giving a flag of -no_ex will prevent using extra rotamer \
        sampling during repacking, saving a lot of time with reduced sampling.')
    parser.add_argument('-no_models', '--make_models', action='store_false', 
        help='Giving a flag of -no_models will prevent PDB models from being \
        generated. This will prevent the energetic calculations and not yield \
        athe substituted models, but will skip the most time-consuming steps.')
    parser.add_argument('-part', '--parallel_partition', type=int, nargs=2, 
        default=[1,1], help='To parallelize the analysis, enter two integers. \
        The first is the number of partitions, the second is which partition \
        member to run on this processor, from 1 to the number of partitions. \
        not working for now.')
    parser.add_argument('-ind', '--ind_type', choices = ['pose', 'pdb'], \
        default = 'pdb', help='To show mutation residues indices in pdb or in pose order')
    parser.add_argument('-op', '--only_protein', action='store_true', 
        help='Giving a flag of -op will prevent ligands and RNA motifs from repacking.')
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


def compare_sequences(pdb_name, pdb_seq, seq_2, query_seq, ind_by): #ind_by has two option, pose or pdb
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

    hyphen_list2 = find(seq_2_for_model, "-")
    if len(hyphen_list2) != 0:
        for xx in hyphen_list2:
            seq_2_for_model[xx] = seq_1_for_model[xx]

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
    Input a 3-member tuple representing a point substitution, of the form: 
    (site, original_AA, substituted_AA)
    Note: this is the output of the compare_sequences function.

    Returns the typical string summary of a point substitution, such as S328A, 
    which indicates that site 328 mutated from S to A.
    """
    m =     point_substitution[4] + \
            point_substitution[1] + \
            str(point_substitution[0]) + \
            point_substitution[2]

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
    check_pose = pose.split_by_chain()[main_chain]
    
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


def check_coevolution(pose, score_function, substitutions):
    """
    Given a pose and a score function, and a list of substitutions, makes 
    selectors for the substitution sites and checkes whether the selections
    have nonzero interaction energy.

    substitutions should be list of 3-member tuples, each representing a point 
    substitution in the form output by the compare_sequences function: 
    (site, original_AA, substituted_AA)    
    """
    # Make a list of substituted sites
    substitution_sites = [s[0] for s in substitutions]

    # Set output to False unless an interaction is found.
    sites_interact = False

    # Automatically no interaction if there are not multiple substitutions
    if len(substitution_sites) < 2:
        return sites_interact

    # Non-redundantly iterate pairwise through the list of substitutions 
    # Skip last res to avoid self-comparison
    for n, sub1 in enumerate(substitution_sites[:-2]):
        # Make selector for first residue 
        res1 = ResidueIndexSelector(str(sub1))

        # Iterate through all subsequent substitution sites
        for sub2 in substitution_sites[n + 1:]:
            # Make selector for second residue 
            res2 = ResidueIndexSelector(str(sub2))

            # Determine interface energy
            intE = interaction_energy(pose, score_function, res1, res2)

            # If the energy is nonzero, residues are interacting
            if intE != 0:
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


def make_point_mutant_task_factory(site_changes, ex12=True, repacking_range=False):
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
# Initialize selection of all altered residues
    # Make TaskFactory to input changes
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    if ex12:
        tf.push_back(ExtraRotamers(0, 1, 1))
        tf.push_back(ExtraRotamers(0, 2, 1))

    # Force packing to new residue for each mutated site
    mutated_res_selection = OrResidueSelector() # Keep list of mobile residues
    for pm in site_changes:
        res_selection = ResidueIndexSelector(str(pm[0]))
        aa_force = RestrictAbsentCanonicalAASRLT()
        aa_force.aas_to_keep(pm[2])
        tf.push_back(OperateOnResidueSubset(aa_force, res_selection))
        mutated_res_selection.add_residue_selector(res_selection)
    if repacking_range:
        # Repack surrounding residues within [repacking_range] angstroms
        repacking_res_selection = NeighborhoodResidueSelector()
        repacking_res_selection.set_focus_selector(mutated_res_selection)
        repacking_res_selection.set_distance(repacking_range)
        repacking_res_selection.set_include_focus_in_subset(False)
        repack = RestrictToRepackingRLT()
        tf.push_back(OperateOnResidueSubset(repack, repacking_res_selection))
        # Do nothing to other residues
        mutated_and_repacking_res_selection = OrResidueSelector(mutated_res_selection, repacking_res_selection)
        prevent = PreventRepackingRLT()
        # By setting the arg2 to True, the mutation and repacking res selection is flipped
        tf.push_back(OperateOnResidueSubset(prevent, mutated_and_repacking_res_selection, True))
    else:
        # Repack all other residues to accommodate substitutions
        repack = RestrictToRepackingRLT()
        # By setting the arg2 to True, the mutation selection is flipped
        tf.push_back(OperateOnResidueSubset(repack, mutated_res_selection, True))

    return tf


def make_point_mutant_task_factory_only_protein(site_changes, ex12=True, \
    repacking_range=False):
    """
    Only repacking protein residues.
    Written by Zhuofan
    """
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    if ex12:
        tf.push_back(ExtraRotamers(0, 1, 1))
        tf.push_back(ExtraRotamers(0, 2, 1))

    mutated_res_selection = OrResidueSelector()
    for pm in site_changes:
        res_selection = ResidueIndexSelector(str(pm[0]))
        aa_force = RestrictAbsentCanonicalAASRLT()
        aa_force.aas_to_keep(pm[2])
        tf.push_back(OperateOnResidueSubset(aa_force, res_selection))
        mutated_res_selection.add_residue_selector(res_selection)

    protein_selection = ResiduePropertySelector(ResidueProperty.PROTEIN)

    if repacking_range:
        neighboring_res_selection = NeighborhoodResidueSelector()
        neighboring_res_selection.set_focus_selector(mutated_res_selection)
        neighboring_res_selection.set_distance(repacking_range)
        neighboring_res_selection.set_include_focus_in_subset(False)

        repacking_res_selection = AndResidueSelector(\
            neighboring_res_selection, protein_selection)

        repack = RestrictToRepackingRLT()
        tf.push_back(OperateOnResidueSubset(repack, repacking_res_selection))

        mutated_and_repacking_res_selection = OrResidueSelector(\
            mutated_res_selection, repacking_res_selection)
        prevent = PreventRepackingRLT()
        tf.push_back(OperateOnResidueSubset(prevent, mutated_and_repacking_res_selection, \
            True))

    else:
        not_portein_selection = NotResidueSelector(protein_selection)
        mutated_and_static_res_selection = OrResidueSelector(mutated_res_selection, \
            not_portein_selection)

        repack = RestrictToRepackingRLT()
        tf.push_back(OperateOnResidueSubset(repack, mutated_and_static_res_selection, True))

        prevent = PreventRepackingRLT()
        tf.push_back(OperateOnResidueSubset(prevent, not_portein_selection))

    return tf


def make_move_map(pose, only_protein=False):
    """
    Generates a movemap for either (Minimization after Repacking) or 
    (only FastRelax protein scaffold).
    Written by Zhuofan
    """
    move_map = pr.MoveMap()
    move_map.set_bb(True)
    if only_protein:
        protein_selection = ResiduePropertySelector(ResidueProperty.PROTEIN)
        protein_res_true_vector = protein_selection.apply(pose)
        move_map.set_chi(protein_res_true_vector)
    else:
        move_map.set_chi(True)
    move_map.set_jump(True)
    return move_map


def repacking_with_muts_and_minimization(pose, task_factory, move_map, score_function):
    """
    Applies point mutations to a given pose. This is done through a 
    PackRotamersMover, followed by minimization.
    pose is a Rosetta pose
    site_changes is a list of point_mutation objects
    Modified by Zhuofan.
    """
    # Make a copy Pose
    mutated_pose = pr.Pose(pose)

    # Repacking
    if task_factory:
        # Apply changes with PackRotamersMover
        prm = PackRotamersMover()
        prm.score_function(score_function)
        prm.task_factory(task_factory)

        # Apply the PackRotamersMover
        prm.apply(mutated_pose)

    # Minimization after repacking
    # If there is no mutation (i.e., the reference sequence), just run minimization
    min_mover = MinMover()
    min_mover.movemap(move_map)
    min_mover.score_function(score_function)

    # Apply the MinMover to the modified Pose
    min_mover.apply(mutated_pose)

    return mutated_pose


def fast_relax_with_muts(pose, task_factory, move_map, score_function, decoys):
    """
    Applies point mutations to a given pose. This is done through a 
    Fast Relax Mover.
    Written by Zhuofan
    """
    # Make FastRelax mover
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    fast_relax.set_movemap(move_map)

    traj_idx = 0
    if_first_decoy = True
    while traj_idx < decoys:
        pose_copy = pr.Pose()
        pose_copy.assign(pose)
        fast_relax.apply(pose_copy)
        decoy_energy = total_energy(pose_copy, score_function)
        if if_first_decoy:
            if_first_decoy = False
            mutated_pose = pose_copy
            lowest_energy = decoy_energy
        elif decoy_energy < lowest_energy:
            mutated_pose = pose_copy
            lowest_energy = decoy_energy
        traj_idx += 1

    return mutated_pose


def make_mutant_model(ref_pose, substitutions, score_function, ex12=True, \
    repacking_range=False, relax_decoys=False, only_protein=False, \
    debugging_mode=False):
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
    # Make a task factory to substitute residues and repack
    if len(substitutions) > 0:
        if only_protein:
            tf = make_point_mutant_task_factory_only_protein(substitutions, \
                ex12=ex12, repacking_range=repacking_range)
        else:
            tf = make_point_mutant_task_factory(substitutions, ex12=ex12, \
                repacking_range=repacking_range)
        if debugging_mode:
            print(tf.create_task_and_apply_taskoperations(ref_pose))
    else:
        tf = None

    # Set up a flexible-backbone movemap
    mm = make_move_map(ref_pose, only_protein=only_protein)

    # Make residue changes
    if len(substitutions) == 0 or (not relax_decoys):
        # If there is no mutation (i.e., the reference sequence), just 
        # run a minimization. Otherwise run both repacking and minimization.
        mutated_pose = repacking_with_muts_and_minimization(ref_pose, \
            tf, mm, score_function)
    else:
        mutated_pose = fast_relax_with_muts(ref_pose, tf, mm, \
            score_function, relax_decoys)

    # Initialize data collection dict
    mutated_pose_data = {}

    # Add energy change to output data
    wt_energy = total_energy(ref_pose, score_function)
    mutant_energy = total_energy(mutated_pose, score_function)
    mutated_pose_data['energy_change'] = mutant_energy - wt_energy

    # Add RMSD to output data
    mutated_pose_data['rmsd'] = get_rmsd(ref_pose, mutated_pose)

    return mutated_pose, mutated_pose_data

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
    ['DKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIV \
    QLSEISMDNSPNLAWPLIVTALRA',840], ['KMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQG', 954]]
    """
    wild_seqs = {}
    tmp_pose = wild_pose.clone()
    remove_nonprotein_residues(tmp_pose)
    if len(list_fasta_names) > 1:
        cut = cut.strip().split(",")
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
 
def analyze_mutant_protein(seqrecord, ref_pose, sf, query, pdb_seq, fa_ind, pdb_name, ind_type, main_chain=1,  
    make_model=True, ex12=True, oligo_chains=None, cat_res=None, substrate_chains=None,
    cut_order = None, rep_fa_ind = None, replicate_id1 = None,
    repacking_range=False, relax_decoys=False, only_protein=False, debugging_mode=False):
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
    substitutions, new_subs = compare_sequences(pdb_name, {list(fa_ind.keys())[0] : pdb_seq[list(fa_ind.keys())[0]]}, seqrecord.seq, 
        query[list(fa_ind.values())[0]], ind_type)

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
                    rest.seq, query[fa_inds[ff]], ind_type)
               substitutions += res_substitutions
               new_subs += res_new_subs
            replicate_id1.pop(mut_tags['id_1']) 
            
    # Add substitutions to output data
    sub_shorts = [get_mutant_brief(pm) for pm in substitutions]
    mut_tags['substitutions'] = ';'.join(sub_shorts)

    # Add whether mutation is in pose or not
    mut_tags['is_in_PDB'] = ';'.join([str(pm[3]) for pm in substitutions])

    # Add to output data whether there are multiple substitutions
    mut_tags['multiple'] = len(new_subs) > 1

    # Add whether substitution(s) are conservative to output data
    conservations = [is_conservative(pm) for pm in substitutions]
    mut_tags['conservative'] = ';'.join([str(i) for i in conservations])

    # If the mutant model is being made, the subsequent steps should use the 
    # Mutant model for increased accuracy, so that step is performed here
    if make_model:
        mutated_pose, mut_pose_info = \
            make_mutant_model(ref_pose, new_subs, sf,  ex12=ex12, \
                repacking_range=repacking_range, relax_decoys=relax_decoys, \
                only_protein=only_protein, debugging_mode=debugging_mode)

        # Add mutant model info to the output data
        mut_tags.update(mut_pose_info)

        # Check if is_in_pdb, replace NA with energy terms
        tmp = mut_tags['is_in_PDB'].split(";")
        if 'True' not in tmp:
            mut_tags['energy_change'] = 'NA'
            mut_tags['rmsd'] = 'NA'

        # Set future calculations to use the mutated pose
        pose = mutated_pose
    else:
        mutated_pose = None
        pose = ref_pose

    # Add substitution layer to output data
    layers = [identify_res_layer(pose, pm[0], main_chain=main_chain) 
        for pm in new_subs]
    mut_tags['layer'] = ';'.join(layers)

    # Add whether multiple substitutions are interacting to output data
    mut_tags['coevolution'] = check_coevolution(pose, sf, new_subs)

    # Add whether substitutions are near multi-chain interface to output data
    if oligo_chains:
        interface_checks = [check_interface_proximity(pose, sf, pm[0], 
            sym_chains, main_chain=main_chain) for pm in new_subs]
        mut_tags['near_interface'] = ';'.join([str(i) for i in interface_checks])

    # Add whether substitutions are near catalytic residues to output data
    if cat_res:
        cat_checks = [check_catalytic_proximity(pose, sf, pm[0], cat_res)
            for pm in new_subs]
        mut_tags['near_catalytic'] = ';'.join([str(i) for i in cat_checks])

    # Add whether substitutions are near substrate to output data
    if substrate_chains:
        substrate_checks = [check_interface_proximity(pose, sf, pm[0], 
            substrate_chains, main_chain=main_chain) for pm in new_subs]
        mut_tags['near_substrate'] = ';'.join([str(i) for i in substrate_checks])

    # Convert to a DataFrame
    mut_df = pd.DataFrame(mut_tags, index=[1])

    return mut_df, mutated_pose, substitutions, new_subs, replicate_id1

########## Executing ###########################################################         

def main(args):
    # Load Rosetta
    opts = '-mute all -run:preserve_header'
    if args.params:
        opts += ' -extra_res_fa ' + ' '.join(args.params)
    pr.init(opts)

    # Set up output directory if generating PDB models
    outdir = out_directory(args.out_dir)

    # Load protease model
    wild_pose = pr.pose_from_pdb(args.template_pdb)

    #Set membrane scorefunction if needed
    if args.membrane:
        base_sf = 'franklin2019'
    else:
        base_sf = 'ref2015_cst'

    # Get score function, applying symmetry if specified
    if args.symmetry:
        sf = SymmetricScoreFunction()
        sf.add_weights_from_file(base_sf)

        sfsm = SetupForSymmetryMover(args.symmetry)
        sfsm.apply(wild_pose)
    else:
        sf = pr.create_score_function(base_sf)
    
    #Set up membrane for membrane protein
    if args.membrane:
        if args.symmetry:
            add_memb = pr.rosetta.protocols.membrane.symmetry.SymmetricAddMembraneMover(args.span_file)
        else:
            add_memb = pr.rosetta.protocols.membrane.AddMembraneMover(args.span_file)
        add_memb.apply(wild_pose)

    # Apply constraints
    if args.constrain:
        wild_pose = coord_constrain_pose(wild_pose)

    # Read fasta sequences
    # edited by Changpeng
    # Find replicates of reference proteins in multiple FASTA files
    # And read wild_seq for each monomer

    list_fasta_names = args.mutants_list.strip().split(",")
    print("------------------------------")
    print("Considered FASTA files are:")
    print(list_fasta_names)

    replicates = check_replicates(list_fasta_names)

    #Read all FASTA files. need to combine with check_replicates together further
    analyze_lists = []
    wts = []
    for i in list_fasta_names:
        fasta_list = parse_fastafile(i)
        analyze_lists.append(fasta_list[1:])
        wts.append(fasta_list[0])
    # get the replicate inds in FASTA files
    replicate_inds = replicate_seqs(replicates,analyze_lists)
 
    # Take subset of fasta list if parallelizing
    #analyze_list = partition_list(fasta_list[1:], *args.parallel_partition)
    #print(analyze_list[0])
        
    #wt = fasta_list[0]
    # split wild_pose into different parts that are corresponding to cut_regions
    pdb_seqs = cut_by_chain(wild_pose, args.cut_region_by_chains, list_fasta_names)
    print("Cut regions of the pdb and their start index in the pdb are:")
    print(pdb_seqs)

    if args.cut_region_by_chains is not None:
        cut = args.cut_region_by_chains.strip().split(",")
    else:
        cut = list(pdb_seqs.keys())[0]

    # Iterate through all identified fasta sequences, altering protease model
    all_mutants_info = pd.DataFrame([])
    all_substitutions = []
    for c in range(len(list_fasta_names)):
        for n, mutant in enumerate(analyze_lists[c]):
            single_mutant_info, mutated_pose, substitutions, new_subs, replicate_inds = \
            analyze_mutant_protein(mutant, wild_pose, sf, wts, pdb_seqs, {cut[c]: c},
                pdb_name = args.template_pdb,  
                main_chain=args.main_chain, make_model=args.make_models, 
                ex12=args.extra_rotamers, oligo_chains=args.interface_chain, 
                cat_res=args.catalytic_residues, 
                substrate_chains=args.ligand_chain,
                repacking_range=args.neighborhood_residue,
                relax_decoys=args.fast_relax, debugging_mode=args.debugging_mode,
                ind_type = args.ind_type, cut_order = cut, rep_fa_ind = replicates, replicate_id1 = replicate_inds, 
                only_protein=args.only_protein)

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
                    if sub[3] == True:
                        subs_clean.append((sub[4],sub[0],sub[1],sub[2]))
                mut_summary = [get_mutant_brief(s) for s in substitutions]
                str_mut_sum = '-'.join(mut_summary)
                outname = '{}_{}_{}.pdb'.format(tag, first_id, str_mut_sum)
                outname = join(outdir, outname)
                mutated_pose.dump_pdb(outname)

            # Append individual mutant data to aggregate
            all_mutants_info = all_mutants_info.append(single_mutant_info)
            for s in subs_clean:
                if s not in all_substitutions:
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
    all_subs_info = pd.DataFrame(all_substitutions)
    all_subs_info.columns = ['chain', 'site', 'native', 'mutant']
    all_subs_info = all_subs_info.sort_values(by=['chain', 'site', 'mutant'])
    all_subs_info.to_csv(subs_report_name, index=False)


if __name__ == '__main__':
    args = parse_args()
    main(args)
