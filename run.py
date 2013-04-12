# This example file shows how to melt an SiO2 system

from dsmcconfig import *

# n{x,y,z} is number of processors in each dimension.
# unit_cells_{x,y,z} is number of unit cells in each dimension.
# More parameters in constructor.

program = DSMC()
dsmc = program.compile(skip_compile=False)

program.reset()
program.prepare_new_system()

program.run_dsmc()
program.reservoir_fraction = 0

program.timesteps = 10000
program.pressure_A = 200000
program.pressure_B = 100000
program.create_config_file()

program.run_dsmc()
