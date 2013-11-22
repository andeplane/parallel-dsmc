from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.005, nx=1, ny=1, nz=1)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.compile(skip_compile=False, opengl_support=True)
geometry.save_file = True
geometry.visualize = False
geometry.num_voxels_x = 128
geometry.num_voxels_y = 128
geometry.num_voxels_z = 128
geometry.binary_output_folder = "../worlds/sinus/"
program.flow_direction = 2
program.atoms_per_molecule = 100
program.Lx = 1.0
program.Ly = 1.0
program.Lz = 1.0
program.seed = 10011126

program.temperature = 300

program.world = geometry.binary_output_folder
program.velocity_bins = 128

amplitude = 0.025
displacement = 0.4
actual_Ly = program.Ly - 2*(amplitude + displacement)
program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=actual_Ly)
program.cells_x = 16
program.cells_y = 16
program.cells_z = 16

for sinus_mode in [0,1,2,4,8]:
    program.reset()
    
    #geometry.create_sinus(sinus_mode=4, amplitude=amplitude, displacement=displacement)    
    program.prepare_new_system()
    program.run(dsmc)
    program.statistics_interval = 1000
    program.timesteps = 20000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="states/00_thermalized")

if True:
    program.statistics_interval = 5
    program.timesteps = 100000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="states/01_sampling")
