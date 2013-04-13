# mdconfig.py
# Module to simplify MD simulations. 
# If logging is enabled, every terminal command is written to run_log.txt
# 
# Anders Hafreager, 2013

import subprocess
import os
import logging
from datetime import datetime

class DSMC:
	def __init__(self, compiler = "icpc", dt=0.001, logging_enabled=True):
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
		self.constants = dict()
		self.constants['boltzmann'] = 1.3806488e-23

		self.load_previous_state = False
		self.create_movie = False

		self.atoms_per_molecule = 1000
		self.timesteps = 10000
		self.temperature = 300
		self.wall_temperature = 300

		self.statistics_interval = 100
		self.movie_every_n_frame = 1
		self.movie_molecules = 10000

		# Gravity direction also determins which dimension we measure flux in
		self.gravity = 0
		self.gravity_direction = 2

		self.maintain_pressure = True
		self.pressure_A = 100000
		self.pressure_B = 100000
		self.reservoir_fraction = 0

		self.world = "../worlds/empty.bin"

		# physical constants, density in N/micrometer^3, length in micrometer, mass in amu, viscosity in Pa*s
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
		self.logging_enabled = logging_enabled
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
		self.run_command("cp -r "+path+" ./")

	def save_state(self, path):
		if self.test_mode: return
		self.log("Saving state to "+str(path))
		self.run_command("mkdir -p "+path)
		self.run_command("cp -r log state_files statistics run_log.txt dsmc.ini Tocontinue "+path)

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
		self.run_command('rm -rf state_files')
		self.run_command('rm dsmc.ini')
		self.run_command('rm Tocontinue')
		
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
		self.run_command('echo 0.0 0 > Tocontinue')

	def compile(self, path = "./program/dsmc", skip_compile = False):
		if not skip_compile:
			current_directory = os.getcwd()
			os.chdir(path)
			self.run_command('qmake')
			self.run_command('make clean')
			self.run_command('make')
			os.chdir(current_directory)

			self.run_command('mv '+path+'/dsmc ./dsmc')
			

		if not os.path.isfile("./dsmc"):
			print "Executable ./dsmc is not compiled, aborting!"
		return './dsmc'

	def create_config_file(self, config_file='dsmc_original.ini'):
		"""
		Creates config file with current settings.
		
		"""
		
		self.log("Creating config file")
		original_file = open("program/000_config_files/"+config_file,'r');
		output_file   = open('dsmc.ini','w');

		for line in original_file:
			line = line.replace('__load_previous_state__',str(self.load_previous_state).lower() )
			line = line.replace('__create_movie__',str(self.create_movie).lower() )
			line = line.replace('__atoms_per_molecule__',str(self.atoms_per_molecule) )
			line = line.replace('__timesteps__',str(self.timesteps) )
			line = line.replace('__temperature__',str(self.temperature) )
			line = line.replace('__wall_temperature__',str(self.wall_temperature) )
			line = line.replace('__dt__',str(self.dt) )
			line = line.replace('__statistics_interval__',str(self.statistics_interval) )
			line = line.replace('__movie_every_n_frame__',str(self.movie_every_n_frame) )
			line = line.replace('__movie_molecules__',str(self.movie_molecules) )
			line = line.replace('__gravity__',str(self.gravity) )
			line = line.replace('__gravity_direction__',str(self.gravity_direction) )
			line = line.replace('__maintain_pressure__',str(self.maintain_pressure).lower() )
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
	
	def run_dsmc(self):
		"""
		Runs specified executable, puts data into a folder, named <project_name>-<name>.
		
		:param executable: Specifies which executable to run.
		:type executable: str.
		:param name: Specifies which output data to use. Default is the last data output.
		:type name: str.
		
		"""

		self.total_timesteps += self.timesteps
		if self.test_mode: return
		executable = "dsmc"
		if not os.path.isfile(executable):
			print "Executable "+executable+" is not compiled, aborting!"
			exit()
		
		self.log("Running executable "+executable)
		now = datetime.now()
		self.run_command("./"+executable+" | tee log")
		t1 = (datetime.now() - now).seconds
		steps_per_second = self.timesteps / max(t1,1)

		self.log("Process used %d seconds (%d timesteps per second)" % ( t1, steps_per_second ))

	def prepare_new_system(self):
		self.timesteps = 1
		self.load_previous_state = False
		self.create_config_file()
		self.load_previous_state = True
