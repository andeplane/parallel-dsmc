from dsmcconfig import *
from dsmc_unit_converter import *
from dsmcplotting import *
from math import sqrt

program = DSMC(dt=0.001, nx=8, ny=8, nz=8)
dsmc = program.compile(skip_compile=True, name="job1")
uc = DSMC_unit_converter(program)

program.reset()

program.gravity = 0.01
program.flow_direction = 2
program.atoms_per_molecule = 10
program.Lx = 1.0
program.Ly = 1.0
program.Lz = 1.0
program.seed = 100

actual_ly = program.Ly*0.9

program.density = uc.density_from_knudsen_number(knudsen_number=0.01, length=actual_ly)

program.temperature = 300
program.cells_x = 120
program.cells_y = 120
program.cells_z = 120

program.world = "../worlds/box_0.9_201.bin"
program.velocity_bins = 200
program.velocity_profile_type = "box"

if True:
    program.prepare_new_system()
    program.run(dsmc)
    program.statistics_interval = 1000
    program.timesteps = 10000
    program.create_config_file()
    program.run(dsmc)

    program.save_state(path="states/00_thermalized")

if True:
    program.load_state(path="states/00_thermalized")

    program.statistics_interval = 100
    program.timesteps = 100000
    program.create_config_file()
    program.run(dsmc)
    program.save_state(path="states/01_sampling")

program.mkdir('plots')
vis = Visualizer()
#vis.plot_velocity_distribution_box(state_path="./", height=program.Ly, output = "plots/velocity_profile.pdf", show_plot=False)
x,v,p = vis.polyfit_velocity("states/01_sampling", height=1.0)

mean_free_path=0.009e-6
a = p[0]
b = p[1]
c = p[2]
dvdx = 2*a*x[0] + b
alpha = v[0]/(mean_free_path*dvdx)
print "dvdx(x_slip) = ", dvdx 
print "alpha= ", alpha