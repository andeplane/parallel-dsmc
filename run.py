# This example file shows how to melt an SiO2 system

from dsmcconfig import *

# n{x,y,z} is number of processors in each dimension.
# unit_cells_{x,y,z} is number of unit cells in each dimension.
# More parameters in constructor.

program = DSMC()
dsmc = program.compile(skip_compile=False, name="job1")

program.reset()

program.reservoir_fraction = 0.2
program.atoms_per_molecule = 100

#program.density = 1e25
program.temperature = 1000

program.Lx  = 0.1
program.Ly  = 0.1
program.Lz  = 0.1
program.nodes_x = 1
program.nodes_y = 1
program.nodes_z = 1
program.cells_per_node_x = 10
program.cells_per_node_y = 10
program.cells_per_node_z = 10

ideal_gas_pressure = program.density*program.constants['boltzmann']*program.temperature;
program.pressure_A = ideal_gas_pressure + 200000
program.pressure_B = ideal_gas_pressure
program.world = "../worlds/box.bin"

program.prepare_new_system()
program.run(dsmc)

program.create_xyz("./state_files/", xyz_file = "./state.xyz")

#program.timesteps = 1000
#program.create_config_file()
#program.run(dsmc)
