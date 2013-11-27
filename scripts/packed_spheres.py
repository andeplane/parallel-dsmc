from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.005, nx=1, ny=1, nz=1)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.binary_output_folder = "../worlds/packed_spheres/"
program.mkdir(geometry.binary_output_folder)

program.atoms_per_molecule = 1

program.seed = 1001

program.world = geometry.binary_output_folder
program.velocity_bins = 128
program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=program.Ly)

for spheres_num in range(40,500,20):
	program.reset()

	state_base_folder = "states/%d" % (spheres_num)
	program.mkdir(state_base_folder)
	geometry.create_packed_spheres(radius = 0.05, spheres_num = spheres_num, inverted=True)

	num_particles = program.get_number_of_particles(geometry)
	#if(num_particles > 10000000): program.atoms_per_molecule = program.atoms_per_molecule*10

	program.apply_pressure_gradient_percentage(factor=1.1)
	program.set_number_of_cells(geometry, particles_per_cell=10)

	program.prepare_new_system()
	program.run(dsmc)
	program.statistics_interval = 10
	program.timesteps = 20000
	program.create_config_file()
	program.run(dsmc)
	program.save_state(path="%s/00_thermalized" % (state_base_folder))

	program.statistics_interval = 5
	program.timesteps = 100000
	program.create_config_file()
	program.run(dsmc)
	program.save_state(path="%s/01_sampling" % (state_base_folder))
	program.run_command("cp %s/porosity.txt %s" % (geometry.binary_output_folder, state_base_folder))