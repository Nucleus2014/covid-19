import argparse
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task.operation import \
	ExtraRotamers, IncludeCurrent, RestrictToRepacking
from pyrosetta.rosetta.core.scoring import ScoreType as st
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from os import makedirs
from os.path import basename, isdir, join

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
		help="What do you want to name the relaxed PDB? (Default appends \
		'_relaxed' to original name.)")
	parser.add_argument('-n', "--n_decoys", type=int, default=100, 
		help="How many decoys do you want? (Default: 100)")
	parser.add_argument('-sf', "--score_function", default='ref2015_cst',
		help="What score function do you want to use? (Default: ref2015_cst)")
	parser.add_argument('-cat', "--cat_res", type=int, default=None, nargs='+',
		help="List residues that should be immobile, separated by spaces")
	parser.add_argument('-cst', "--constraints", default=None,
		help="If EnzDes constraints are to be applied in addition to the \
		default coordinate constraints, specify the file")
	parser.add_argument('-lig', "--ligand", default=None,
		help="If there is a ligand, specify the params file")
	parser.add_argument('-sym', "--symmetry", default=None,
		help="If the relax should be symmetric, specify the symdef file")
	parser.add_argument('-ccw', "--coord_wt", type=float, default=None,
		help="Specify the coordinate constraints weight (Default: 1.0)")
	parser.add_argument('-edw', "--enzdes_wt", type=float, default=None,
		help="Specify the constraints weight for enzdes constraints \
		(Default: 1.0)")
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

	# Creating symmetry setup
	if args.symmetry:
		sfsm = SetupForSymmetryMover(args.symmetry)

	# Setting up the scorefunction with the desired constraint weights
	if args.symmetry:
		sf = SymmetricScoreFunction()
		sf.add_weights_from_file(args.score_function)
	else:
		sf = create_score_function(args.score_function)

	if args.coord_wt:
		sf.set_weight(st.coordinate_constraint, args.coord_wt)

	if args.enzdes_wt:
		sf.set_weight(st.atom_pair_constraint, args.enzdes_wt)
		sf.set_weight(st.angle_constraint, args.enzdes_wt)
		sf.set_weight(st.dihedral_constraint, args.enzdes_wt)

	# Creating FastRelax protocol with the given score function
	fr = FastRelax()
	fr.set_scorefxn(sf)

	# Packer tasks
	tf = standard_task_factory()
	tf.push_back(RestrictToRepacking())
	tf.push_back(IncludeCurrent())
	tf.push_back(ExtraRotamers(0, 1, 1))
	tf.push_back(ExtraRotamers(0, 2, 1))
	fr.set_task_factory(tf)

	# Determining file name
	if args.name: 
		file_name = args.name
	else:
		file_name = basename(args.pdb_file).replace('.pdb', '').replace('.gz', '')
		file_name += '_relaxed'

	out_name = join(dir_name, file_name)

	# Loading PDB file, applying constraints
	pose = pose_from_pdb(args.pdb_file)
	ac.apply(pose)
	if args.constraints:
		enz_cst.apply(pose)
	if args.symmetry:
		sfsm.apply(pose)

	# RMSD metric
	rmsdm = RMSDMetric()
	rmsdm.set_comparison_pose(pose)

	# Running relax set
	jd = PyJobDistributor(out_name, args.n_decoys, sf)
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
