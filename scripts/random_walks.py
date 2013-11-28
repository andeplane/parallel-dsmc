from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.001, nx=2, ny=2, nz=2)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

random_walk_run_id = "1"

geometry = DSMC_geometry(program)
geometry.binary_output_folder = "../worlds/random_walk_%s/" % (random_walk_run_id)
program.mkdir(geometry.binary_output_folder)
program.atoms_per_molecule = 1
program.Lx = 0.5
program.Ly = 0.5
program.Lz = 0.5

program.seed = 1001

program.world = geometry.binary_output_folder
program.velocity_bins = 128
program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=program.Ly)


for seed in range(1,10,1):
	program.reset()

	state_base_folder = "../random_walk_states/%s-%d" % (random_walk_run_id, seed)
	program.mkdir(state_base_folder)
	geometry.visualize = False
	geometry.save_state = True
	geometry.create_random_walk(seed=seed, inverted = False, walker_number = 100, walker_steps = 300, walker_turn_probability = 1, walker_thickness_change_prob = 0, walker_max_thickness = 5)
	
	num_particles = program.get_number_of_particles(geometry)
	
	program.apply_pressure_gradient_percentage(factor=1.1)
	program.set_number_of_cells(geometry, particles_per_cell=10)

	program.prepare_new_system()
	program.run(dsmc)
	program.statistics_interval = 1000
	program.timesteps = 20000
	program.create_config_file()
	program.run(dsmc)
	program.save_state(path="%s/00_thermalized" % (state_base_folder))
	
	program.statistics_interval = 1000
	program.timesteps = 100000
	program.create_config_file()
	program.run(dsmc)
	program.save_state(path="%s/01_sampling" % (state_base_folder))
	program.run_command("cp %s/porosity.txt %s" % (geometry.binary_output_folder, state_base_folder))