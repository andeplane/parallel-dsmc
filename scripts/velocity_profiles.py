from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.01, nx=2, ny=2, nz=2)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.compile(skip_compile=False, opengl_support=True)
geometry.save_file = True
geometry.visualize = False
geometry.num_voxels_x = 128
geometry.num_voxels_y = 128
geometry.num_voxels_z = 128
geometry.binary_output_folder = "../worlds/poiseuille/"
geometry.create_poiseuille()

program.atoms_per_molecule = 1
program.Lx = 1.0
program.Ly = 1.0
program.Lz = 1.0
program.seed = 1001

k = 10.0
kn = 2.0/sqrt(pi)*k

program.density = uc.density_from_knudsen_number(knudsen_number=kn, length=program.Ly)
program.temperature = 300
program.apply_pressure_gradient_percentage(factor=1.1)
program.set_number_of_cells(geometry, particles_per_cell=10)

program.world = geometry.binary_output_folder
program.velocity_bins = 100

if True:
    program.reset()
    program.prepare_new_system()
    program.run(dsmc)
    program.statistics_interval = 1000
    program.timesteps = 20000
    program.create_config_file()
    program.run(dsmc)

    program.save_state(path="states/00_thermalized")

if True:
    program.statistics_interval = 5
    program.timesteps = 200000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="states/01_sampling")
