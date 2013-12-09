from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.005, nx=2, ny=2, nz=2)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.binary_output_folder = "../worlds/one_cylinder/"
program.mkdir(geometry.binary_output_folder)
program.atoms_per_molecule = 100
program.world = geometry.binary_output_folder
geometry.create_cylinders(radius=0.45, num_cylinders_per_dimension = 1)

knudsen_numbers = [0.01]+[i*0.1 for i in range(101)]
for kn in knudsen_numbers:
    state_folder = "states/%.6f/" % (kn)
    program.reset()
    program.density = uc.density_from_knudsen_number(knudsen_number=0.01, length=program.Ly)
    num_particles = program.get_number_of_particles(geometry)
    if(num_particles < 250000 and program.atoms_per_molecule > 1) program.atoms_per_molecule = program.atoms_per_molecule/10
    program.apply_pressure_gradient_percentage(factor=1.1)
    program.set_number_of_cells(geometry, particles_per_cell=20)

    program.prepare_new_system()
    program.run(dsmc)
    program.statistics_interval = 1000
    program.timesteps = 20000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="%s/00_thermalized" % (state_folder))

    program.statistics_interval = 1000
    program.timesteps = 200000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="%s/01_sampling" % (state_folder))