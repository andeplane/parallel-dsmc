from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.005, nx=4, ny=4, nz=4)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.binary_output_folder = "../worlds/packed_spheres/"
program.mkdir(geometry.binary_output_folder)

program.world = geometry.binary_output_folder
program.velocity_bins = 128
program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=program.Ly)

#for spheres_radius in [0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]:
for spheres_radius in [0.01, 0.015, 0.02 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.065, 0.07, 0.075, 0.08]:
	program.reset()

	state_base_folder = "states/%f" % (spheres_radius)
	program.mkdir(state_base_folder)
	geometry.create_packed_spheres(radius = spheres_radius, spheres_type = 1, wanted_porosity = 0.3)

	num_particles = program.get_number_of_particles(geometry)
	program.atoms_per_molecule = 10

	program.apply_pressure_gradient_percentage(factor=1.1)
	program.set_number_of_cells(geometry, particles_per_cell=20)

	program.prepare_new_system()
	program.run(dsmc)
	program.statistics_interval = 10
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