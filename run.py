# This example file shows how to melt an SiO2 system

from dsmcconfig import *

# n{x,y,z} is number of processors in each dimension.
# unit_cells_{x,y,z} is number of unit cells in each dimension.
# More parameters in constructor.

program = DSMC()
dsmc = program.compile(skip_compile=False, name="job1")

program.reset()
program.reservoir_fraction = 0.2
program.atoms_per_molecule = 500

#program.density = 1e25
program.wall_temperature = 100
program.temperature = 100
program.cells_x = 20
program.cells_y = 20
program.cells_z = 20

ideal_gas_pressure = program.density*program.constants['boltzmann']*program.temperature;
program.maintain_pressure = True
program.pressure_A = ideal_gas_pressure + 50000
program.pressure_B = ideal_gas_pressure
program.world = "../worlds/box.bin"
#program.surface_interaction = "cercignani_lampis"

program.prepare_new_system()
program.run(dsmc)

program.timesteps = 20000
program.create_config_file()
program.run(dsmc)

program.run(dsmc)