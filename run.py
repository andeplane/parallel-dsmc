# This example file shows how to melt an SiO2 system

from dsmcconfig import *

# n{x,y,z} is number of processors in each dimension.
# unit_cells_{x,y,z} is number of unit cells in each dimension.
# More parameters in constructor.

program = DSMC()
dsmc = program.compile(skip_compile=False)

for density in [1e25, 1e26, 1e27]:
    for atoms_per_molecule in [10000,5000,2000,1000,500,250,100]:
        program.reset()

        program.density = density
        program.temperature = 300
        ideal_gas_pressure = program.density*program.constants['boltzmann']*program.temperature;

        program.world = "../worlds/box.bin"
        program.prepare_new_system()

        program.run_dsmc()
        program.reservoir_fraction = 0.6
        program.atoms_per_molecule = atoms_per_molecule
        program.timesteps = 10000
        program.pressure_A = ideal_gas_pressure + 100000
        program.pressure_B = ideal_gas_pressure
        program.create_config_file()

        program.run_dsmc()

        program.timesteps = 250000
        program.create_config_file()

        program.run_dsmc()
        program.save_state("states/%e-%d" % (density,atoms_per_molecule))