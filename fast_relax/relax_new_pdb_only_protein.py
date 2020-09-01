import argparse
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.pack.task.operation import \
    IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
    RestrictAbsentCanonicalAASRLT, RestrictToRepackingRLT, \
    PreventRepackingRLT
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.select.residue_selector import \
    ResiduePropertySelector
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.membrane import AddMembraneMover
from pyrosetta.rosetta.protocols.membrane.symmetry import \
    SymmetricAddMembraneMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from os import makedirs
from os.path import basename, isdir, join
import sys

'''
When downloading a new PDB file, relax with coordinate constraints to eliminate clashes.
Requires a PDB file input.
Options:
Name (-n, string): change the output PDB name from [original_name]_relaxed.pdb
Score function (-sf, string): change the score function from the default of ref2015_cst
Catalytic residues (-cat, int, multiple accepted): list residues that should not be moved 
'''

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file", help="What PDB file do you want to relax?")
    parser.add_argument("-od", "--out_dir", 
        help="Name an output directory for decoys (Default: current directory)")
    parser.add_argument('-name', "--name", 
        help="What do you want to name the relaxed PDB? (Default appends '_relaxed' to original name.)")
    parser.add_argument('-n', "--n_decoys", type=int, default=100, 
        help="How many decoys do you want? (Default: 100)")
    parser.add_argument('-rep', "--repulsive_type", choices=['hard', 'soft'],
        default='hard', help="What score function do you want to use?")
    parser.add_argument('-symm', "--symmetry", type=str, 
        help='Symmetric file.')
    parser.add_argument('-memb', '--membrane', required=False, action='store_true',
        help='Declare if the protein is a membrane protein.')
    parser.add_argument('-mspan', '--span_file', 
        required=any(x in ['--membrane','-memb'] for x in sys.argv),
        help='If the pose is a membrane protein, include a spanfile.')
    parser.add_argument('-cat', "--cat_res", type=int, default=None, nargs='+',
        help="List residues that should be immobile, separated by spaces")
    parser.add_argument('-cst', "--constraints", default=None,
        help="If EnzDes constraints are to be applied in addition to the \
        default coordinate constraints, specify the file")
    parser.add_argument('-lig', "--ligand", default=None,
        help="If there is a ligand, specify the params file")
    parser.add_argument('-ccw', "--coord_wt", type=float, default=1,
        help="Specify the coordinate constraints weight (Default: 1.0)")
    parser.add_argument('-edw', "--enzdes_wt", type=float, default=None,
        help="Specify the constraints weight for enzdes constraints (Default: 1.0)")
    args = parser.parse_args()
    return args


def main(args):
    # Destination folder for PDB files
    if args.out_dir:
        dir_name = args.out_dir
        if not isdir(dir_name):
            makedirs(dir_name)
    else:
        dir_name = ""

    # Creating coordinate constraints for the entire molecule
    cg = CoordinateConstraintGenerator()
    ac = AddConstraints()
    ac.add_generator(cg)

    # Create enzdes constraints
    if args.constraints:
        enz_cst = AddOrRemoveMatchCsts()
        enz_cst.set_cst_action(ADD_NEW)

    ''' Declare the score function. '''
    if args.symmetry: # Declare symmetric score functions
        score_function = SymmetricScoreFunction()
        if args.repulsive_type == 'hard':
            if args.membrane:
                score_function.add_weights_from_file('franklin2019')
            else:
                score_function.add_weights_from_file('ref2015')
        elif args.repulsive_type == 'soft':
            if args.membrane: # Set up a soft-rep version of franklin2019 manually
                score_function.add_weights_from_file('ref2015_soft')
                score_function.set_weight(ScoreType.fa_water_to_bilayer, 1.0)
            else:
                score_function.add_weights_from_file('ref2015_soft')
    else: # Declare ordinary score functions
        if args.repulsive_type == 'hard':
            if args.membrane:
                score_function = create_score_function('franklin2019')
            else:
                score_function = create_score_function('ref2015')
        elif args.repulsive_type == 'soft':
            if args.membrane: # Set up a soft-rep version of franklin2019 manually
                score_function = create_score_function('ref2015_soft')
                score_function.set_weight(ScoreType.fa_water_to_bilayer, 1.0)
            else:
                score_function = create_score_function('ref2015_soft')

    if args.coord_wt:
        score_function.set_weight(ScoreType.coordinate_constraint, args.coord_wt)

    if args.enzdes_wt:
        score_function.set_weight(ScoreType.atom_pair_constraint, args.enzdes_wt)
        score_function.set_weight(ScoreType.angle_constraint, args.enzdes_wt)
        score_function.set_weight(ScoreType.dihedral_constraint, args.enzdes_wt)

    # Loading PDB file
    pose = pose_from_pdb(args.pdb_file)

    if args.symmetry: # Applying symmetry if specified
        sfsm = SetupForSymmetryMover(args.symmetry)
        sfsm.apply(pose)
        if args.membrane: # Set up membrane for membrane protein
            add_memb = SymmetricAddMembraneMover(args.span_file)
            add_memb.apply(pose)
    else:
        if args.membrane: # Set up membrane for membrane protein
            add_memb = AddMembraneMover(args.span_file)
            add_memb.apply(pose)

    # Creating FastRelax protocol with the given score function
    fr = FastRelax()
    fr.set_scorefxn(score_function)

    # Packer tasks
    tf = standard_task_factory()
    tf.push_back(IncludeCurrent())
    tf.push_back(ExtraRotamers(0, 1, 1))
    tf.push_back(ExtraRotamers(0, 2, 1))
    protein_selector = ResiduePropertySelector(ResidueProperty.PROTEIN)
    repack = RestrictToRepackingRLT()
    tf.push_back(OperateOnResidueSubset(repack, protein_selector))
    prevent = PreventRepackingRLT()
    tf.push_back(OperateOnResidueSubset(prevent, protein_selector, True))
    fr.set_task_factory(tf)

    move_map = MoveMap()
    if args.repulsive_type == 'hard':
        move_map.set_bb(True)
    elif args.repulsive_type == 'soft':
        ''' When using the soft-rep score function, backbone should be fixed. '''
        move_map.set_bb(False)
    protein_res_true_vector = protein_selector.apply(pose)
    move_map.set_chi(protein_res_true_vector)
    fr.set_movemap(move_map)

    # Determining file name
    if args.name: 
        file_name = args.name
    else:
        file_name = basename(args.pdb_file).replace('.pdb', '_relaxed')

    out_name = join(dir_name, file_name)

    # Applying constraints
    
    ac.apply(pose)
    if args.constraints:
        enz_cst.apply(pose)

    # RMSD metric
    rmsdm = RMSDMetric()
    rmsdm.set_comparison_pose(pose)

    print(tf.create_task_and_apply_taskoperations(pose))

    # Running relax set
    jd = PyJobDistributor(out_name, args.n_decoys, score_function)
    while not jd.job_complete:
        pp = Pose()
        pp.assign(pose)
        fr.apply(pp)
        rmsdm.apply(pp)
        jd.output_decoy(pp)


if __name__ == '__main__':
    args = parse_args()

    opts = '-ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false -cst_fa_weight 1.0'
    if args.constraints:
        opts += ' -enzdes::cstfile {} -run:preserve_header'.format(args.constraints)
    if args.ligand:
        opts += ' -extra_res_fa {}'.format(args.ligand)
    init(opts)
    
    main(args)

# '-relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains'
