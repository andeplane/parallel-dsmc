# mdconfig.py
# Module to simplify MD simulations. 
# If logging is enabled, every terminal command is written to run_log.txt
# 
# Anders Hafreager, 2013

import subprocess
import os
import logging
from datetime import datetime
from math import sqrt
from dsmc_unit_converter import *

class DSMC:
	def __init__(self, compiler = "icpc", dt=0.001, nx=1, ny=1, nz=1):
		"""
		Initializes an object with parameters.
		"""
		current_directory = os.getcwd()
		# if os.path.basename(current_directory) == "base_code":
		# 	# We don't want to modify the base_code folder
		# 	print "Cannot run code from base_code folder. Please copy the program into another folder. Aborting!"
		# 	exit()

		self.dt = dt
		self.compiler = compiler
		self.mpi_hostfile = None
		self.constants = dict()
		self.constants['boltzmann'] = 1.3806488e-23

		self.load_previous_state = False
		self.create_movie_files = False
		self.surface_interaction = "thermal"
		self.alpha_t = 0.8
		self.alpha_n = 0.5

		self.max_molecules_per_node = 1e6
		self.atoms_per_molecule = 1000
		self.timesteps = 10000
		self.temperature = 300
		self.wall_temperature = 300

		self.nx = nx
		self.ny = ny
		self.nz = nz

		self.seed = -1
		self.statistics_interval = 100
		self.sampling_bins = 128
		self.velocity_profile_type = "other"
		self.movie_every_n_frame = 1
		self.movie_molecules = 10000

		# Gravity direction also determins which dimension we measure flux in
		self.gravity = 0
		self.flow_direction = 2

		self.pressure_A = 100000
		self.pressure_B = 100000
		self.reservoir_fraction = 0

		self.world = "../worlds/empty.bin"
		
		# physical constants, density in N/m^3, length in micrometer, mass in amu, viscosity in Pa*s
		self.density = 2.4143e25
		self.diam = 3.62e-4
		self.mass = 39.948
		self.viscosity = 2.23e-5

		# Grid size (measured in micrometer)
		self.Lx  = 1
		self.Ly  = 1
		self.Lz  = 1

		self.cells_x = 10
		self.cells_y = 10
		self.cells_z = 10
		
		self.test_mode = False
		self.logging_enabled = True
		self.total_timesteps = 0

		logging.basicConfig(filename='run_log.txt',level=logging.INFO)
		self.log("System initialized")
		
	def log(self, message):
		if self.logging_enabled:
			logging.info("#LOG "+str(datetime.now())+" : "+message)

	def load_state(self, path):
		if self.test_mode: return
		self.clean()
		self.log("Loading state from "+path)
		self.run_command("cp -r "+path+"/ ./")
		self.load_previous_state = True

	def save_state(self, path):
		if self.test_mode: return
		self.log("Saving state to "+str(path))
		self.run_command("mkdir -p "+path)
		self.run_command("cp -r log state_files statistics run_log.txt dsmc.ini Tocontinue "+path)
		self.run_command("cp run.py "+path+"/run_script.py")

	def run_command(self, cmd):
		"""
		Runs command <cmd> in terminal.
		
		:param cmd: Command to run.
		:type cmd: str.
		
		"""
		if self.logging_enabled:
			logging.info("#CMD "+str( datetime.now() )+" : "+cmd)
		subprocess.call(cmd, shell=True)

	def clean(self):
		if self.test_mode: return

		self.run_command('rm log')
		self.run_command('rm -rf state_files/*')
		self.run_command('rm dsmc.ini')
		self.run_command('rm Tocontinue')

	def mkdir(self, directory_name):
		self.run_command('mkdir '+directory_name)
		
	def reset(self):
		"""
		Deletes files created by the md program.
		
		"""
		if self.test_mode: return

		self.log('Resetting root folder')
		self.clean()
		self.run_command('mkdir state_files')
		self.run_command('mkdir movie_files')
		self.run_command('mkdir statistics')
		self.run_command('mkdir statistics/time')
		self.run_command('echo 0.0 0 0 0 > Tocontinue')

	def compile(self, path = "./program/dsmc", name="dsmc", skip_compile = False, clean = False):
		if not skip_compile:
			current_directory = os.getcwd()
			os.chdir(path)
			self.run_command('qmake')
			if clean: self.run_command('make clean')
			self.run_command('make')
			os.chdir(current_directory)
			move_command = 'mv '+path+'/dsmc ./%s' % name
			self.run_command(move_command)
			

		if not os.path.isfile("./"+name):
			print "Executable ./"+name+" is not compiled, aborting!"
			exit()
		return './%s' % name

	def create_config_file(self, config_file='dsmc_original.ini'):
		"""
		Creates config file with current settings.
		
		"""
		
		self.log("Creating config file")
		original_file = open("program/000_config_files/"+config_file,'r');
		output_file   = open('dsmc.ini','w');

		for line in original_file:
			line = line.replace('__load_previous_state__',str(self.load_previous_state).lower() )
			line = line.replace('__create_movie__',str(self.create_movie_files).lower() )
			line = line.replace('__seed__',str(self.seed) )
			line = line.replace('__nx__',str(self.nx) )
			line = line.replace('__ny__',str(self.ny) )
			line = line.replace('__nz__',str(self.nz) )
			line = line.replace('__max_molecules_per_node__',str(self.max_molecules_per_node) )
			line = line.replace('__atoms_per_molecule__',str(self.atoms_per_molecule) )
			line = line.replace('__timesteps__',str(self.timesteps) )
			line = line.replace('__surface_interaction_model__',str(self.surface_interaction).lower() )
			line = line.replace('__alpha_n__',str(self.alpha_n) )
			line = line.replace('__alpha_t__',str(self.alpha_t) )
			line = line.replace('__sampling_bins__', str(self.sampling_bins))
			line = line.replace('__surface_interaction_model__',str(self.surface_interaction).lower() )
			line = line.replace('__velocity_profile_type__', str(self.velocity_profile_type).lower() )
			line = line.replace('__temperature__',str(self.temperature) )
			line = line.replace('__wall_temperature__',str(self.wall_temperature) )
			line = line.replace('__dt__',str(self.dt) )
			line = line.replace('__statistics_interval__',str(self.statistics_interval) )
			line = line.replace('__movie_every_n_frame__',str(self.movie_every_n_frame) )
			line = line.replace('__movie_molecules__',str(self.movie_molecules) )
			line = line.replace('__gravity__',str(self.gravity) )
			line = line.replace('__flow_direction__',str(self.flow_direction) )
			line = line.replace('__pressure_A__',str(self.pressure_A) )
			line = line.replace('__pressure_B__',str(self.pressure_B) )
			line = line.replace('__reservoir_fraction__',str(self.reservoir_fraction) )
			line = line.replace('__world__', str(self.world) )
			line = line.replace('__density__', str(self.density) )
			line = line.replace('__diam__', str(self.diam) )
			line = line.replace('__mass__', str(self.mass) )
			line = line.replace('__viscosity__', str(self.viscosity) )
			line = line.replace('__Lx__', str(self.Lx) )
			line = line.replace('__Ly__', str(self.Ly) )
			line = line.replace('__Lz__', str(self.Lz) )
			line = line.replace('__cells_x__', str(self.cells_x) )
			line = line.replace('__cells_y__', str(self.cells_y) )
			line = line.replace('__cells_z__', str(self.cells_z) )
			
			output_file.writelines(line)

		original_file.close()
		output_file.close()
	
	def run(self, executable="dsmc"):
		"""
		Runs specified executable, puts data into a folder, named <project_name>-<name>.
		
		:param executable: Specifies which executable to run.
		:type executable: str.
		:param name: Specifies which output data to use. Default is the last data output.
		:type name: str.
		
		"""

		self.total_timesteps += self.timesteps
		if self.test_mode: return
		
		if not os.path.isfile(executable):
			print "Executable "+executable+" is not compiled, aborting!"
			exit()
		
		self.log("Running executable "+executable)
		now = datetime.now()
		num_procs = self.nx*self.ny*self.nz
		if self.mpi_hostfile is None:
			self.run_command("mpirun -n %d %s | tee log" % (num_procs, executable))
		else:
			self.run_command("mpirun -n %d --hostfile %s %s | tee log" % (num_procs, self.mpi_hostfile, executable))
		t1 = (datetime.now() - now).seconds
		steps_per_second = self.timesteps / max(t1,1)

		self.log("Process used %d seconds (%d timesteps per second)" % ( t1, steps_per_second ))

	def prepare_new_system(self):
		self.timesteps = 1
		self.load_previous_state = False
		self.create_config_file()
		self.load_previous_state = True

	def create_movie(self, frames):
		num_procs = self.nx*self.ny*self.nz
		self.run_command("%s -O3 program/tools/create_movie.cpp -o create_movie" % self.compiler)
		self.run_command("./create_movie ./ %d %d" % (num_procs, frames) )

	def create_xyz(self, state = "./", xyz_file = "./state.xyz"):
		if self.test_mode: return
		state = state+"state_files/"
		self.run_command(self.compiler + " program/tools/create_xyz.cpp -o ./create_xyz")
		self.run_command("./create_xyz %d %s %s" % (1, state, xyz_file))
	
	def get_mean_free_path(self):
		pi = 3.141592653586
		return 1.0/(sqrt(2)*pi*self.diam**2*self.density*1e-18)*1e-6

	def validate_number_of_cells(self, num_voxels_x, num_voxels_y, num_voxels_z):
		num_voxels_per_cell_x = float(num_voxels_x) / self.cells_x
		num_voxels_per_cell_y = float(num_voxels_y) / self.cells_y
		num_voxels_per_cell_z = float(num_voxels_z) / self.cells_z
		if(num_voxels_per_cell_x % 1 > 0): return False
		if(num_voxels_per_cell_y % 1 > 0): return False
		if(num_voxels_per_cell_z % 1 > 0): return False

		num_cells_per_cpu_x = float(self.cells_x)/self.nx
		num_cells_per_cpu_y = float(self.cells_y)/self.ny
		num_cells_per_cpu_z = float(self.cells_z)/self.nz

		if(num_cells_per_cpu_x % 1 > 0): return False
		if(num_cells_per_cpu_y % 1 > 0): return False
		if(num_cells_per_cpu_z % 1 > 0): return False

		return True
	def apply_pressure_gradient_percentage(self, factor):
		uc = DSMC_unit_converter(self)
		ideal_gas_pressure = self.density*self.constants["boltzmann"]*self.temperature
		pressure_A = ideal_gas_pressure*factor
		pressure_B = ideal_gas_pressure
		delta_p = pressure_A - pressure_B
		self.gravity = uc.acceleration_from_si(uc.pressure_difference_to_gravity(delta_p, self.Lz))

	def get_number_of_particles(self, geometry):
		uc = DSMC_unit_converter(self)
		volume = self.Lx*self.Ly*self.Lz*geometry.get_porosity()
		# Find number of particles to set correct cell count
		num_particles =  int(uc.number_density_from_si(self.density)*volume / self.atoms_per_molecule)
		return num_particles

	def set_number_of_cells(self, geometry, particles_per_cell=10):
		num_particles = self.get_number_of_particles(geometry)
		print "Estimated number of particles: ", num_particles
		porosity = geometry.get_porosity()
		number_of_cells = num_particles / particles_per_cell * porosity
		number_of_cells_per_dimension = int(number_of_cells**(1.0/3.0)) # Assuming cubi system
		
		self.cells_x = max(number_of_cells_per_dimension,1)
		self.cells_y = max(number_of_cells_per_dimension,1)
		self.cells_z = max(number_of_cells_per_dimension,1)

		self.cells_x = self.cells_x - self.cells_x%self.nx
		self.cells_y = self.cells_y - self.cells_y%self.ny
		self.cells_z = self.cells_z - self.cells_z%self.nz

		while(geometry.num_voxels_x%self.cells_x>0): 
			self.cells_x += self.nx
			self.cells_x = self.cells_x - self.cells_x%self.nx
			if(self.cells_x > geometry.num_voxels_x): 
				print "Number of cells (%d) exceeds number of voxels (%d), choosing cells_x=num_voxels_x=%d!" % (self.cells_x, geometry.num_voxels_x, geometry.num_voxels_x)
				self.cells_x = geometry.num_voxels_x
				break
		while(geometry.num_voxels_y%self.cells_y>0): 
			self.cells_y += self.ny
			self.cells_y = self.cells_y - self.cells_y%self.ny
			if(self.cells_y > geometry.num_voxels_y): 
				print "Number of cells (%d) exceeds number of voxels (%d), choosing cells_y=num_voxels_y=%d!" % (self.cells_y, geometry.num_voxels_y, geometry.num_voxels_y)
				self.cells_y = geometry.num_voxels_y
				break
		while(geometry.num_voxels_z%self.cells_z>0): 
			self.cells_z += self.nz
			self.cells_z = self.cells_z - self.cells_z%self.nz
			if(self.cells_z > geometry.num_voxels_z): 
				print "Number of cells (%d) exceeds number of voxels (%d), choosing cells_x=num_voxels_x=%d!" % (self.cells_z, geometry.num_voxels_z, geometry.num_voxels_z)
				self.cells_z = geometry.num_voxels_z
				break

		print "Cells: (%d, %d, %d)" % (self.cells_x, self.cells_y, self.cells_z)
		num_cells = self.cells_x*self.cells_y*self.cells_z
		particles_per_cell = num_particles / num_cells
		print "Particles per cell: ", particles_per_cell