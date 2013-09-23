from dsmcconfig import *
from dsmcplotting import *
from dsmc_unit_converter import *

from math import sqrt

program = DSMC()
dsmc = program.compile(skip_compile=False, name="job1")
uc = DSMC_unit_converter(program)

program.reset()

program.gravity = 0.0
program.gravity_direction = 2
program.reservoir_fraction = 0.2
program.atoms_per_molecule = 50
program.Lx = 0.2
program.Ly = 0.2
program.Lz = 1.0

actual_ly = program.Ly*0.8 # The rest is upper and lower wall

program.density = uc.density_from_knudsen_number(knudsen_number=0.1, length=actual_ly)

program.temperature = 300
program.cells_x = 12
program.cells_y = 12
program.cells_z = 60

tube_length = program.Lz*(1.0 - program.reservoir_fraction)

pressure_difference = uc.gravity_to_pressure_difference(g=0.01, length=tube_length)
program.maintain_pressure = True
ideal_gas_pressure = program.density*program.constants['boltzmann']*program.temperature;
program.pressure_A = ideal_gas_pressure + pressure_difference
program.pressure_B = ideal_gas_pressure

print "P_A = ", program.pressure_A
print "P_B = ", program.pressure_B

program.world = "../worlds/box_fraction_0.2.bin"
program.statistics_interval = 1000
program.velocity_bins = 80

if True:
	program.prepare_new_system()
	program.run(dsmc)

	program.timesteps = 5000
	program.create_config_file()
	program.run(dsmc)
	program.save_state(path="states/00_thermalized")

	program.load_state(path="states/00_thermalized")
	program.statistics_interval = 10
	program.timesteps = 100000

	program.create_config_file()
	program.run(dsmc)
	program.save_state(path="states/01_thermalized")

if False:
	program.load_state(path="states/01_thermalized")
	program.create_movie_files = True
	program.movie_molecules = 8000
	program.timesteps = 1000
	program.create_config_file()
	program.run(dsmc)
	program.create_movie(frames=1000)
mass_si = 6.63352088e-26; #kg
factor = 1.0#/sqrt(2*program.constants['boltzmann']/mass_si*program.temperature)

vis = Visualizer()
vis.plot_velocity_distribution_box(state_path="states/01_thermalized", height=program.Ly, output="velocity_profile.pdf", show_plot=True, factor=factor)