import os

from dsmcconfig import *

class DSMC_geometry:
	def __init__(self, dsmc):
		self.dsmc = dsmc
		self.type = "box"
		self.screen_width = 1680
		self.screen_height = 1050
		self.num_voxels_x = 128
		self.num_voxels_y = 128
		self.num_voxels_z = 128
		self.number_of_neighbor_averages = 1
		self.save_file = True
		self.binary_output_folder = "../worlds/test/"
		self.num_processors_x = dsmc.nx
		self.num_processors_y = dsmc.ny
		self.num_processors_z = dsmc.nz
		self.create_border = False

		# Sphere settings
		self.sphere_radius = 0.8
		self.inverted = True

		# Cylinder settings
		self.cylinder_radius = 0.8
		self.num_cylinders_per_dimension = 1

		# Diamond square settings
		self.hurst_exponent = 0.5
		self.diamond_square_h0 = 1.0

		# Sinus settings
		self.amplitude = 1.0
		self.sinus_mode = 1
		self.displacement = 1.0

		# Random walk settings
		self.seed = -1
		self.walker_number = 5
		self.walker_steps = 1000
		self.walker_max_thickness = 5
		self.walker_thickness_change_prob = 0.03
		self.walker_turn_probability = 0.1

		# Perlin settings
		self.perlin_octave = 1
		self.perlin_frequency = 1
		self.perlin_amplitude = 1
		self.perlin_num_scales = 4
		self.perlin_constant = 8.143352246
		self.perlin_scale_factor = 3.35156131
		self.perlin_threshold = 0.5

		# Visualizer settings
		self.Lx = 1.0
		self.Ly = 1.0
		self.Lz = 1.0
		self.scale = 1
		self.visualize = False
		self.marching_cubes_threshold = 0.5

	def create_box(self):
		self.type = "box"
		self.create_config_file()
		self.run()

	def create_poiseuille(self):
		self.type = "poiseuille"
		self.create_config_file()
		self.run()

	def create_sphere(self, radius = 0.9, inverse=True):
		self.inverted = inverse
		self.sphere_radius = radius
		self.type = "sphere"
		self.create_config_file()
		self.run()

	def create_cylinders(self, radius = 0.9, num_cylinders_per_dimension = 1):
		self.cylinder_radius = radius
		self.num_cylinders_per_dimension = num_cylinders_per_dimension
		self.type = "cylinders"
		self.create_config_file()
		self.run()

	def create_diamond_square(self, hurst_exponent = 0.5, seed = 1, diamond_square_h0 = 1):
		self.hurst_exponent = hurst_exponent
		self.perlin_seed = seed
		self.diamond_square_h0 = diamond_square_h0
		self.type = "diamond_square"
		self.create_config_file()
		self.run()

	def create_empty(self):
		self.type = "empty"
		self.create_config_file()
		self.run()

	def create_random_walk(self, seed=-1, inverted = False, walker_number = 1, walker_steps = 100, walker_turn_probability = 0.1, walker_thickness_change_prob = 0.01, walker_max_thickness = 5):
		self.type = "random_walk"
		self.seed = seed
		self.walker_turn_probability = walker_turn_probability
		self.walker_number = walker_number
		self.walker_steps = walker_steps
		self.walker_thickness_change_prob = walker_thickness_change_prob
		self.walker_max_thickness = walker_max_thickness
		self.inverted = inverted
		self.create_config_file()
		self.run()

	def create_sinus(self, sinus_mode=1, amplitude=1.0, displacement=1.0):
		self.type = "sinus"
		self.sinus_mode = sinus_mode
		self.amplitude = amplitude
		self.displacement = displacement
		self.create_config_file()
		self.run()

	def create_perlin(self, threshold=0.5, num_scales = 4, scale_factor=3.64341524563, constant=8.23552246134, seed = 1, amplitude=1, frequency=1, octave=1):
		self.perlin_threshold = threshold
		self.perlin_constant = constant
		self.seed = seed
		self.perlin_amplitude = amplitude
		self.perlin_frequency = frequency
		self.perlin_octave = octave
		self.perlin_scale_factor = scale_factor
		self.perlin_num_scales = num_scales
		self.type = "perlin"
		self.create_config_file()
		self.run()

	def create_config_file(self, config_file='geometry_settings_original.ini'):
		original_file = open("program/000_config_files/"+config_file,'r');
		output_file   = open('geometry_settings.ini','w');
		for line in original_file:
			line = line.replace('__type__',str(self.type).lower() )
			line = line.replace('__num_voxels_x__',str(self.num_voxels_x) )
			line = line.replace('__num_voxels_y__',str(self.num_voxels_x) )
			line = line.replace('__num_voxels_z__',str(self.num_voxels_x) )
			line = line.replace('__seed__',str(self.seed) )
			line = line.replace('__walker_number__',str(self.walker_number) )
			line = line.replace('__walker_steps__',str(self.walker_steps) )
			line = line.replace('__walker_thickness_change_prob__',str(self.walker_thickness_change_prob) )
			line = line.replace('__walker_max_thickness__',str(self.walker_max_thickness) )
			line = line.replace('__walker_turn_probability__',str(self.walker_turn_probability) )
			line = line.replace('__number_of_neighbor_averages__',str(self.number_of_neighbor_averages) )
			line = line.replace('__diamond_square_h0__',str(self.diamond_square_h0) )
			line = line.replace('__save_file__',str(self.save_file).lower() )
			line = line.replace('__binary_output_folder__',str(self.binary_output_folder) )
			line = line.replace('__num_processors_x__',str(self.num_processors_x) )
			line = line.replace('__num_processors_y__',str(self.num_processors_y) )
			line = line.replace('__num_processors_z__',str(self.num_processors_z) )
			line = line.replace('__sphere_radius__', str(self.sphere_radius))
			line = line.replace('__amplitude__',str(self.amplitude) )
			line = line.replace('__displacement__',str(self.displacement) )
			line = line.replace('__sinus_mode__', str(self.sinus_mode))
			line = line.replace('__cylinder_radius__', str(self.cylinder_radius))
			line = line.replace('__num_cylinders_per_dimension__', str(self.num_cylinders_per_dimension))
			line = line.replace('__inverted__',str(self.inverted).lower() )
			line = line.replace('__create_border__',str(self.create_border).lower() )
			line = line.replace('__perlin_octave__', str(self.perlin_octave) )
			line = line.replace('__perlin_frequency__',str(self.perlin_frequency) )
			line = line.replace('__perlin_amplitude__',str(self.perlin_amplitude) )
			line = line.replace('__perlin_num_scales__',str(self.perlin_num_scales) )
			line = line.replace('__perlin_constant__',str(self.perlin_constant) )
			line = line.replace('__hurst_exponent__',str(self.hurst_exponent) )
			line = line.replace('__perlin_scale_factor__',str(self.perlin_scale_factor) )
			line = line.replace('__perlin_threshold__',str(self.perlin_threshold) )
			line = line.replace('__screen_width__',str(self.screen_width) )
			line = line.replace('__screen_height__',str(self.screen_height) )
			line = line.replace('__Lx__',str(self.Lx) )
			line = line.replace('__Ly__',str(self.Ly) )
			line = line.replace('__Lz__',str(self.Lz) )
			line = line.replace('__scale__',str(self.scale) )
			line = line.replace('__visualize__',str(self.visualize).lower() )
			line = line.replace('__marching_cubes_threshold__', str(self.marching_cubes_threshold) )
			
			output_file.writelines(line)

		original_file.close()
		output_file.close()

	def run_command(self, cmd):
		subprocess.call(cmd, shell=True)

	def compile(self, opengl_support=False, skip_compile = False, clean = False):
		if not skip_compile:
			if clean: self.run_command('make clean')
			if opengl_support: self.run_command('make -f make_geometry')
			else: self.run_command('make -f make_geometry_noopengl')			

		if not os.path.isfile("./geometry"):
			print "Executable ./geometry is not compiled, aborting!"
			exit()
		return './geometry'

	def run(self):
		if not os.path.isfile("./geometry"):
			print "Executable ./geometry is not compiled, aborting!"
			exit()
		
		self.run_command("./geometry")

	def get_porosity(self):
		file_pointer = open(self.binary_output_folder+"/porosity.txt",'r');
		porosity = file_pointer.readline()
		return float(porosity)