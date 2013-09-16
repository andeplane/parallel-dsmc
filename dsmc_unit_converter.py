from dsmcconfig import *
from math import sqrt, pi

class DSMC_unit_converter:
	def __init__(self, dsmc):
		self.dsmc = dsmc
		self.m0 = 1.66053886e-27  # si
		self.L0 = 1e-6            # si
		self.E0 = 1.65088e-21     # si
		self.E0ev = 1.0318e-2     # eV
		self.kb = 1.3806503e-23   # si

		self.t0 = self.L0*sqrt(self.m0/self.E0)
		self.F0 = self.E0/self.L0
		self.T0 = self.E0/self.kb
		self.P0 = self.m0/(self.t0**2*self.L0)
		self.v0 = self.L0/self.t0
		self.a0 = self.v0/self.t0
		self.visc0 = self.P0*self.t0
		self.diff0 = self.L0**2/self.t0
		self.perm0 = self.L0**2
		self.number_density0 = 1.0/(self.L0**3)

	def pressure_to_si(self, P): 
		return P*self.P0

	def pressure_from_si(self, P): 
		return P/self.P0

	def temperature_to_si(self, T):
		return T*self.T0

	def temperature_from_si(self, T):
		return T/self.T0

	def mass_to_si(self, m):
		return m*self.m0

	def mass_from_si(self, m):
		return m/self.m0

	def length_to_si(self, L):
		return L*self.L0

	def length_from_si(self, L):
		return L/self.L0

	def force_to_si(self, F):
		return F*self.F0

	def force_from_si(self, F):
		return F/self.F0

	def energy_to_si(self, E):
		return E*self.E0

	def energy_from_si(self, E):
		return E/self.E0

	def time_to_si(self, t):
		return t*self.t0

	def time_from_si(self, t):
		return t/self.t0

	def velocity_to_si(self, v):
		return v*self.v0

	def velocity_from_si(self, v):
		return v/self.v0

	def acceleration_to_si(self, a):
		return a*self.a0

	def acceleration_from_si(self, a):
		return a/self.a0

	def viscosity_to_si(self, v):
		return v*self.visc0

	def viscosity_from_si(self, v):
		return v/self.visc0

	def permeability_to_si(self, k):
		return k*self.perm0

	def permeability_from_si(self, k):
		return k/self.perm0

	def diffusion_to_si(self, d):
		return d*self.diff0

	def diffusion_from_si(self, d):
		return d/self.diff0

	def number_density_to_si(self, d):
		return d*self.number_density0

	def number_density_from_si(self, d):
		return d/self.number_density0

	def gravity_to_pressure_difference(self, g, density, length):
		''' Converts gravity driven force to equivalent pressure difference
		'''
		return self.acceleration_to_si(g)*density*self.length_to_si(length)

	def density_from_knudsen_number(self, knudsen_number, length):
		''' Calculates required density to get a specific Knudsen number through
				Kn=mean_free_path/length
		''' 
		return self.number_density_to_si(1.0/(sqrt(2.0)*pi*self.dsmc.diam**2*knudsen_number*length))
	
	def length_from_knudsen_number(self, knudsen_number):
		''' Calculates required density to get a specific Knudsen number through
				Kn=mean_free_path/length
		''' 
		return 1.0/(sqrt(2.0)*pi*self.dsmc.diam**2*knudsen_number*self.dsmc.density)