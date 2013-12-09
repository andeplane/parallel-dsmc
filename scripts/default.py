from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.005, nx=2, ny=2, nz=2)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.binary_output_folder = "../worlds/poiseuille/"
program.mkdir(geometry.binary_output_folder)
geometry.create_poiseuille()

program.atoms_per_molecule = 10

program.world = geometry.binary_output_folder
program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=program.Ly)

program.reset()
program.apply_pressure_gradient_percentage(factor=1.1)
program.set_number_of_cells(geometry, particles_per_cell=20)

program.prepare_new_system()
program.run(dsmc)
program.statistics_interval = 1000
program.timesteps = 10000
program.create_config_file()
program.run(dsmc)
program.save_state(path="states/00_thermalized")

program.statistics_interval = 100
program.timesteps = 50000
program.create_config_file()
program.run(dsmc)
program.save_state(path="staets/01_sampling")