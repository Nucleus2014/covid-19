import pyrosetta as pr
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
    RestrictAbsentCanonicalAASRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector, \
    ChainSelector, LayerSelector, InterGroupInterfaceByVectorSelector, \
    NotResidueSelector, OrResidueSelector, ResidueIndexSelector, \
    TrueResidueSelector 
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    InteractionEnergyMetric, RMSDMetric, TotalEnergyMetric
from pyrosetta.rosetta.protocols.minimization_packing import \
    MinMover, PackRotamersMover
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover

def compare_sequences(seq_1, seq_2):
    """
    Given a reference sequence and a comparison sequence, identify the sites 
    where the comparison sequence differs from the reference. Returns a list of 
    3-member lists. Each trio represents one point substitution. The three items 
    listed in the trio are the site, the starting residue type, and the mutant 
    residue type. Site numbers are 1-indexed (not 0-indexed as Python default) 
    to work with Rosetta. If the sequences are not the same length, it breaks 
    the script.
    """
    # Make sure sequences are the same length
    assert len(seq_1) == len(seq_2)

    # Initialize dict
    site_changes = []

    # Identify altered sites
    for n, i in enumerate(seq_2):
        if i != seq_1[n]:
            site = n + 1
            seq_1_value = seq_1[n]
            seq_2_value = i
            
            site_changes.append([site, seq_1_value, seq_2_value])

    return site_changes

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

def identify_res_layer(pose, res_number):
    """
    Determines whether a given residue in a pose is in the core, boundary, or 
    surface layer of the protein.
    """
    # If the PDB has multiple chains
    check_pose = pose.split_by_chain()[1]
    
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

def interaction_energy(pose, score_function, selection_1, selection_2):
    """
    Given a Rosetta pose, score function, and two residue selectors, 
    calculates the interaction energy between the selections.
    """
    interact_metric = InteractionEnergyMetric()
    interact_metric.set_scorefunction(score_function)
    interact_metric.set_residue_selectors(selection_1, selection_2)
    return interact_metric.calculate(pose)

def make_point_mutant_task_factory(site_1, change_1, site_2=None, change_2=None, fast=True):
    """
    Given a site integer and a 1-letter residue name (or two sites and two 
    residue names), generates a TaskFactory to repack all residues 
    in a pose with the altered sites repacking to the new sequence.
    """
    # Initialize selection of all altered residues
    modified_residues = OrResidueSelector() # Keep list of mobile residues

    # Make TaskFactory to input changes
    task_factory = TaskFactory()
    task_factory.push_back(IncludeCurrent())
    if not fast:
        task_factory.push_back(ExtraRotamers(0, 1, 1))
        task_factory.push_back(ExtraRotamers(0, 2, 1))

    # Force packing to new residue for the mutated site
    res_selection_1 = ResidueIndexSelector(str(site_1))
    aa_force_1 = RestrictAbsentCanonicalAASRLT()
    aa_force_1.aas_to_keep(change_1)
    task_factory.push_back(OperateOnResidueSubset(aa_force_1, res_selection_1))
    modified_residues.add_residue_selector(res_selection_1)
    
    # Add second site if necessary
    if site_2 and change_2:
        res_selection_2 = ResidueIndexSelector(str(site_2))
        aa_force_2 = RestrictAbsentCanonicalAASRLT()
        aa_force_2.aas_to_keep(change_2)
        task_factory.push_back(OperateOnResidueSubset(aa_force_2, res_selection_2))
        modified_residues.add_residue_selector(res_selection_2)
            
    # Repack all other residues to accommodate substitutions
    unchanged_residues = NotResidueSelector(modified_residues)
    repack = RestrictToRepackingRLT()
    task_factory.push_back(OperateOnResidueSubset(repack, unchanged_residues))

    return task_factory

def make_point_changes(pose, task_factory, score_function):
    """
    Applies point mutations to a given pose. This is done through a 
    PackRotamersMover, followed by minimization.
    Inputs are a Pose, a TaskFactory, and a ScoreFunction
    """
    # Make PackRotamersMover
    pack_rotamers = PackRotamersMover()
    pack_rotamers.score_function(score_function)
    pack_rotamers.task_factory(task_factory)
    
    # Make a copy Pose and apply the PackRotamersMover
    mutated_pose = pr.Pose(pose)
    pack_rotamers.apply(mutated_pose)

    # Set up fixed-backbone movemap for minimization
    movemap = pr.MoveMap()
    movemap.set_bb(True)
    movemap.set_chi(True)
    movemap.set_jump(True)

    # Create the MinMover
    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(score_function)
    
    # Apply the MinMover to the modified Pose
    min_mover.apply(mutated_pose)

    return mutated_pose

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


