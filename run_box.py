from dsmcconfig import *
from dsmcplotting import *
from dsmc_unit_converter import *

from math import sqrt

vis = Visualizer()

program = DSMC(dt=0.001, nx=2, ny=2, nz=1)
dsmc = program.compile(skip_compile=False, name="job1")
uc = DSMC_unit_converter(program)

program.reset()

program.gravity = 0
program.flow_direction = 2
program.atoms_per_molecule = 10
program.Lx = 1
program.Ly = 1
program.Lz = 1
program.seed = 100

actual_ly = program.Ly

program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=actual_ly)

program.temperature = 300
program.cells_x = 20
program.cells_y = 20
program.cells_z = 20

# Pressure calculations
pressure_difference = uc.gravity_to_pressure_difference(g=0, length=program.Lz)
program.maintain_pressure = False
ideal_gas_pressure = program.density*program.constants['boltzmann']*program.temperature;
program.pressure_A = ideal_gas_pressure
program.pressure_B = ideal_gas_pressure

program.world = "../worlds/empty.bin"
#program.world = "../worlds/la_palma.bin"
#program.world = "../worlds/box_with_rotation_1_101_theta_pi_half.bin"
program.statistics_interval = 100
program.velocity_bins = 100
program.velocity_profile_type = "box"

if False:
	program.prepare_new_system()
	program.run(dsmc)

if True:
	program.prepare_new_system()
	program.run(dsmc)

	print "\n\nPressure in reservoir A: ", program.pressure_A
	print "Pressure in reservoir B: ", program.pressure_B
	print "Average pressure: ", (program.pressure_A + program.pressure_B) / 2.0
	print "Density: ", program.density
	print "\n\n"

	program.statistics_interval = 100
	program.timesteps = 20000
	program.create_config_file()
	program.run(dsmc)
	#vis.plot_linear_density_profile(state_path="./", length=program.Lz, output="plots/density_profile.pdf", show_plot=True)
	exit()

	program.save_state(path="states/00_thermalized")
	
	program.load_state(path="states/00_thermalized")
	program.statistics_interval = 10
	program.timesteps = 10000

	program.create_config_file()
	program.run(dsmc)
	program.save_state(path="states/01_thermalized")

if False:
	program.load_state(path="states/01_thermalized")
	program.statistics_interval = 3
	program.timesteps = 50000

	program.create_config_file()
	program.run(dsmc)
	program.save_state(path="states/02_thermalized")

if False:
	#program.load_state(path="states/01_thermalized")
	program.create_movie_files = True
	program.movie_molecules = 1717
	program.timesteps = 1000
	program.movie_every_n_frame = 1
	program.create_config_file()
	program.run(dsmc)
	program.create_movie(frames=1000)
mass_si = 6.63352088e-26; #kg
factor = 1.0#/sqrt(2*program.constants['boltzmann']/mass_si*program.temperature)

show_plot = False
#vis.plot_linear_density_profile(state_path="./", length=program.Lz, output="plots/density_profile.pdf", show_plot=show_plot)
# vis.plot_linear_pressure_profile(state_path="./", length=program.Lz, output="plots/pressure_profile.pdf", show_plot=show_plot)
# vis.plot_velocity_distribution_box(state_path="./", height=program.Ly, output="plots/velocity_profile.pdf", show_plot=show_plot, factor=factor)
# vis.plot_permeability(state_path="./", output="plots/permeability.pdf", show_plot=show_plot)