from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.005, nx=2, ny=2, nz=2)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.compile(skip_compile=False, opengl_support=True)
geometry.save_file = True
geometry.visualize = False
geometry.num_voxels_x = 128
geometry.num_voxels_y = 128
geometry.num_voxels_z = 128
geometry.binary_output_folder = "../worlds/diamond/"
program.flow_direction = 2
program.atoms_per_molecule = 1
program.Lx = 1.0
program.Ly = 1.0
program.Lz = 1.0
program.seed = 10011126

program.temperature = 300

program.world = geometry.binary_output_folder
program.velocity_bins = 128

program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=program.Ly)

for distance in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    state_base_folder = "states/%f" % (distance)
    program.mkdir(state_base_folder)

    program.reset()
    #geometry.create_diamond_square(hurst_exponent = 0.8, seed = 1, diamond_square_distance = distance)
    num_particles = program.get_number_of_particles(geometry)
    if(num_particles > 1000000): program.atoms_per_molecule = program.atoms_per_molecule*10

    program.apply_pressure_gradient_percentage(factor=1.1)
    program.set_number_of_cells(geometry, particles_per_cell=10)

    program.prepare_new_system()
    program.run(dsmc)
    program.statistics_interval = 1000
    program.timesteps = 20000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="%s/00_thermalized" % (state_base_folder))

    program.statistics_interval = 5
    program.timesteps = 100000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="%s/01_sampling" % (state_base_folder))
