from dsmcconfig import *
from dsmc_unit_converter import *
from math import sqrt, pi
from dsmc_geometry_config import *

program = DSMC(dt=0.005, nx=2, ny=2, nz=2)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.binary_output_folder = "../worlds/one_cylinder_varying_radius/"
program.mkdir(geometry.binary_output_folder)
program.atoms_per_molecule = 100
program.world = geometry.binary_output_folder

radii = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
for radius in radii:
    geometry.create_cylinders(radius=radius, num_cylinders_per_dimension = 1)
    state_folder = "states/%f/" % (radius)
    program.reset()
    program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=radius)
    num_particles = program.get_number_of_particles(geometry)
    if(num_particles < 250000 and program.atoms_per_molecule > 1): program.atoms_per_molecule = program.atoms_per_molecule/10
    
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