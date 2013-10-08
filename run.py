from dsmcconfig import *
from dsmc_unit_converter import *
from dsmcplotting import *

program = DSMC(dt=0.001, nx=1, ny=1, nz=1)
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

program.density = uc.density_from_knudsen_number(knudsen_number=1.0, length=program.Ly)
program.temperature = 300
program.cells_x = 30
program.cells_y = 30
program.cells_z = 30

program.world = "../worlds/empty.bin"
program.velocity_bins = 100
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



program.mkdir('plots')
vis = Visualizer()
vis.plot_velocity_distribution_box(state_path="./", height=program.Ly, output = "plots/velocity_profile.pdf", show_plot=True)