from dsmcconfig import *
from dsmc_unit_converter import *

from math import sqrt

program = DSMC(dt=0.001, nx=2, ny=2, nz=1)
dsmc = program.compile(skip_compile=False, name="job1")
uc = DSMC_unit_converter(program)

program.reset()

program.gravity = 0.1
program.flow_direction = 2
program.atoms_per_molecule = 100
program.Lx = 1.0
program.Ly = 1.0
program.Lz = 1.0
program.seed = 100

actual_ly = program.Ly*0.8

program.density = uc.density_from_knudsen_number(knudsen_number=0.01, length=actual_ly)

program.temperature = 300
program.cells_x = 60
program.cells_y = 60
program.cells_z = 60

# Pressure calculations
program.maintain_pressure = False

program.world = "../worlds/box_with_rotation_1_101_theta_pi_half.bin"
program.velocity_bins = 200
program.velocity_profile_type = "box"

if False:
    program.prepare_new_system()
    program.run(dsmc)
    program.statistics_interval = 1000
    program.timesteps = 10000
    program.create_config_file()
    program.run(dsmc)

    program.save_state(path="states/00_thermalized")

if False:
    program.load_state(path="states/00_thermalized")

    program.statistics_interval = 10
    program.timesteps = 20000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="states/01_sampling")

from dsmcplotting import *

program.mkdir('plots')
vis = Visualizer()
vis.plot_velocity_distribution_box(state_path="./", height=program.Ly, output = "plots/velocity_profile.pdf", show_plot=True)
x,v, p = vis.polyfit_velocity("states/01_sampling", height=1.0)

mean_free_path=0.008e-6
a = p[0]
b = p[1]
c = p[2]
dvdx = 2*a*x[0] + b
alpha = v[0]/(mean_free_path*dvdx)
print "dvdx(x_slip) = ", dvdx 
print "alpha= ", alpha