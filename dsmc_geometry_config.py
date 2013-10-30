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
		self.save_file = False
		self.binary_output_folder = "../worlds/test/"
		self.num_processors_x = dsmc.nx
		self.num_processors_y = dsmc.ny
		self.num_processors_z = dsmc.nz

		# Box settings
		self.box_porosity = 0.9

		# Sphere settings
		self.sphere_radius = 0.8
		self.sphere_inverse = True

		# Diamond square settings
		self.hurst_exponent = 0.5
		self.diamond_square_h0 = 1.0

		# Perlin settings
		self.perlin_octave = 1
		self.perlin_frequency = 1
		self.perlin_amplitude = 1
		self.perlin_seed = 1
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

	def create_box(self, porosity):
		self.box_porosity = porosity
		self.type = "box"
		self.create_config_file()
		self.run()
	def create_poiseuille(self):
		self.type = "poiseuille"
		self.create_config_file()
		self.run()

	def create_sphere(self, radius = 0.9, inverse=True):
		self.sphere_inverse = inverse
		self.sphere_radius = radius
		self.type = "sphere"
		self.create_config_file()
		self.run()

	def create_diamond_square(self, hurst_exponent = 0.5, seed = 1):
		self.hurst_exponent = hurst_exponent
		self.perlin_seed = seed
		self.type = "diamond_square"
		self.create_config_file()
		self.run()

	def create_empty(self):
		self.type = "empty"
		self.create_config_file()
		self.run()

	def create_perlin(self, threshold=0.5, num_scales = 4, scale_factor=3.64341524563, constant=8.23552246134, seed = 1, amplitude=1, frequency=1, octave=1):
		self.perlin_threshold = threshold
		self.perlin_constant = constant
		self.perlin_seed = seed
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
			line = line.replace('__number_of_neighbor_averages__',str(self.number_of_neighbor_averages) )
			line = line.replace('__diamond_square_h0__',str(self.diamond_square_h0) )
			line = line.replace('__save_file__',str(self.save_file).lower() )
			line = line.replace('__binary_output_folder__',str(self.binary_output_folder) )
			line = line.replace('__num_processors_x__',str(self.num_processors_x) )
			line = line.replace('__num_processors_y__',str(self.num_processors_y) )
			line = line.replace('__num_processors_z__',str(self.num_processors_z) )
			line = line.replace('__box_porosity__',str(self.box_porosity) )
			line = line.replace('__sphere_radius__', str(self.sphere_radius))
			line = line.replace('__sphere_inverse__',str(self.sphere_inverse).lower() )
			line = line.replace('__perlin_octave__', str(self.perlin_octave) )
			line = line.replace('__perlin_frequency__',str(self.perlin_frequency) )
			line = line.replace('__perlin_amplitude__',str(self.perlin_amplitude) )
			line = line.replace('__perlin_seed__',str(self.perlin_seed) )
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
