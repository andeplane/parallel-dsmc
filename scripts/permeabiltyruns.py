from dsmcconfig import *
from dsmc_unit_converter import *
from dsmc_geometry_config import *
from math import sqrt

program = DSMC(dt=0.001, nx=2, ny=2, nz=2)
dsmc = program.compile(skip_compile=True, name="main")
uc = DSMC_unit_converter(program)

geometry = DSMC_geometry(program)
geometry.compile(skip_compile=False, opengl_support=False)
geometry.save_file = True
geometry.visualize = False
geometry.num_voxels_x = 128
geometry.num_voxels_y = 128
geometry.num_voxels_z = 128
geometry.binary_output_folder = "../worlds/poiseuille/"
#geometry.create_poiseuille()
porosity = 0.984375

program.flow_direction = 2
program.atoms_per_molecule = 1
program.Lx = 1.0
program.Ly = 1.0
program.Lz = 1.0
program.seed = 1000

program.max_molecules_per_node = 2.5e6
program.temperature = 300

# Pressure calculations
program.maintain_pressure = False
program.world = geometry.binary_output_folder
program.velocity_bins = 128
program.velocity_profile_type = "box"

ideal_gas_pressure = program.density*program.constants["boltzmann"]*program.temperature
pressure_A = ideal_gas_pressure*1.1
pressure_B = ideal_gas_pressure
delta_p = pressure_A - pressure_B
program.gravity = uc.acceleration_from_si(uc.pressure_difference_to_gravity(delta_p, program.Lz))
print "Gravity: ", program.gravity

volume = program.Lx*program.Ly*program.Lz*porosity

for kn in [0.1]:
    state_base_name = 'states/kn=%f' % (kn)
    program.mkdir(directory_name=state_base_name)

    program.density = uc.density_from_knudsen_number(knudsen_number=kn, length=program.Ly)
    num_particles =  int(uc.number_density_from_si(program.density)*volume / program.atoms_per_molecule)
    number_of_cells = num_particles / 20.0 # ~20 particles per cell
    number_of_cells_per_dimension = int(number_of_cells**(1.0/3.0))
    program.cells_x = number_of_cells_per_dimension
    program.cells_y = number_of_cells_per_dimension
    program.cells_z = number_of_cells_per_dimension
    
    program.cells_x = program.cells_x - program.cells_x%program.nx
    program.cells_y = program.cells_y - program.cells_y%program.ny
    program.cells_z = program.cells_z - program.cells_z%program.nz

    while(geometry.num_voxels_x%program.cells_x>0): program.cells_x += program.nx
    while(geometry.num_voxels_y%program.cells_y>0): program.cells_y += program.ny
    while(geometry.num_voxels_z%program.cells_z>0): program.cells_z += program.nz

    print "Running simulation with Kn=%f" % (kn)
    print "num_particles=%d" % (num_particles)
    print "Number of cells: = (%d,%d,%d)" % (program.cells_x, program.cells_y, program.cells_z)
    program.prepare_new_system()
    program.run(dsmc)
    print "New system was created"
    program.statistics_interval = 1000
    program.timesteps = 5000
    program.create_config_file()
    program.run(dsmc)

    program.save_state(path="%s/00_thermalized" % (state_base_name))
    #program.load_state(path="%s/00_thermalized" % (state_base_name))

    program.statistics_interval = 100
    program.timesteps = 5000000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="%s/01_sampling" % (state_base_name))

from dsmcplotting import *
vis = Visualizer()

vis.plot_velocity_distribution_box(state_path="./", height=program.Ly, output = "plots/velocity_profile.pdf", show_plot=True)