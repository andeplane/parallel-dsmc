import subprocess

subprocess.call('icpc -fast -I. create_world_perlin.cpp -o create_world_perlin', shell = True)
subprocess.call('./create_world_perlin 50 50 50', shell = True)
subprocess.call('./create_world ./perlin_m.bin ./perlin.bin -1 0', shell = True)
subprocess.call('../../../dsmc_geometry_visualizer/main', shell = True)