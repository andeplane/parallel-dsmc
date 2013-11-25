from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *
from dsmcplotting import *

program = DSMC(dt=0.001, nx=2, ny=2, nz=2)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.compile(skip_compile=False, opengl_support=True)
geometry.save_file = True
geometry.num_voxels_x = 128
geometry.num_voxels_y = 128
geometry.num_voxels_z = 128
geometry.binary_output_folder = "../worlds/poiseuille/"
geometry.create_poiseuille()

program.atoms_per_molecule = 50
program.Lx = 1.0
program.Ly = 1.0
program.Lz = 1.0
program.seed = 10011126

program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=program.Ly)
program.temperature = 300
program.apply_pressure_gradient_percentage(factor=1.1)

porosity = geometry.get_porosity()
program.set_number_of_cells(geometry, particles_per_cell=10)

geometry.binary_output_folder
program.world = geometry.binary_output_folder
program.sampling_bins = 64

if True:
    program.reset()
    program.prepare_new_system()
    program.run(dsmc)
    program.statistics_interval = 1000
    program.timesteps = 5000
    program.create_config_file()
    program.run(dsmc)

    program.save_state(path="states/00_thermalized")

if True:
    program.statistics_interval = 3
    program.timesteps = 100000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="states/01_sampling")

vis = Visualizer()
vis.plot_velocity_distribution_box(state_path="./", show_plot=True)