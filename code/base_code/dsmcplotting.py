from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from numpy import linspace
from numpy import polyfit
class Visualizer():
	def __init__(self):
		pass

	def plot(self, x, y, x_axis_label = None, y_axis_label = None, curve_label = None, fig = None, filename=None, title=None):
		"""
		Plot y versus x. curve_label will appear in the legend. 
		For multiple curves in one figure, a fig object must be provided on the second and so on call.
		"""
		import matplotlib.pyplot as plt
		if fig is None:
			fig = plt.figure()
	
		plt.figure(fig.number)
		plt.plot(x, y, label=curve_label)
		plt.hold('on')
		if x_axis_label != None:
			plt.xlabel(x_axis_label)
		if y_axis_label != None:	
			plt.ylabel(y_axis_label)
		if curve_label != None:	
			plt.legend()
		if title != None:
			plt.title(title)	
		if filename != None:
			plt.savefig(filename)
		
		return fig

	def combine_two_figures_horizontally(fig1, fig2, combined_fig = None):
		if combined_fig == None:
			combined_fig = "combined.pdf"
		os.system("pdftk %s %s output tmp.pdf" % (fig1, fig2)) 
		os.system("pdfnup --nup 2x1 tmp.pdf" )
		os.system("pdfcrop tmp-nup.pdf ")
		os.system("cp tmp-nup-crop.pdf %s" % combined_fig) 
		os.system("rm tmp.pdf tmp-nup.pdf tmp-nup-crop.pdf")
		return combined_fig

	def plot_velocity_distribution_box(self, state_path="./", height=1.0, output = None, show_plot=False, skip_zeros = True, factor=1.0):
		"""
		Reads a file containing several timesteps with average velocity profile across two parallel plates.
		Function will take average of all timesteps in order to get good statistics.
		"""

		velocities = pylab.loadtxt(state_path+"/statistics/velocity.txt")
		num_bins = len(velocities)   # Number of elements on first line equals number of bins
		velocities = factor*velocities
		print "Plotting velocity profile with "+str(num_bins)+" bins."

		x_axis = linspace(0,height,num_bins) # Create x_axis with possible real channel height
		if skip_zeros:
			self.plot(x=x_axis[velocities>0], y=velocities[velocities>0], x_axis_label="x", y_axis_label="v(x)", filename=output)
		else:
			self.plot(x=x_axis, y=velocities, x_axis_label="x", y_axis_label="v(x)", filename=output)
		if show_plot: plt.show()

	def polyfit_velocity(self, state_path="./", height=1.0):
		"""
		Reads a file containing several timesteps with average velocity profile across two parallel plates.
		Function will take average of all timesteps in order to get good statistics.
		"""

		velocities = pylab.loadtxt(state_path+"/statistics/velocity.txt")
		num_timesteps = len(velocities) # Number of lines equals number of timesteps
		num_bins = len(velocities[0])   # Number of elements on first line equals number of bins
		velocities = sum(velocities,0) / num_timesteps  # Sum each column and normalize to take time average
		
		x_axis = linspace(0,height*1e-6,num_bins) # Create x_axis with possible real channel height
		x_axis = x_axis[velocities>0]
		velocities = velocities[velocities>0]
		p = polyfit(x_axis, velocities, 2)
		return x_axis, velocities, p
		
		
	def plot_linear_density_profile(self, state_path="./", length=1.0, output = None, show_plot=False):
		"""
		Reads a file containing several timesteps with average velocity profile across two parallel plates.
		Function will take average of all timesteps in order to get good statistics.
		"""

		densities = pylab.loadtxt(state_path+"/statistics/linear_density.txt")
		num_timesteps = len(densities) # Number of lines equals number of timesteps
		num_bins = len(densities[0])   # Number of elements on first line equals number of bins
		densities = sum(densities,0) / num_timesteps  # Sum each column and normalize to take time average
		print "Plotting linear density profile with "+str(num_timesteps)+" timesteps and "+str(num_bins)+" bins."

		x_axis = linspace(0,length,num_bins) # Create x_axis with possible real channel height
		
		self.plot(x=x_axis, y=densities, x_axis_label="x", y_axis_label="# molecules", filename=output)
		if show_plot: plt.show()

	def plot_linear_temperature_profile(self, state_path="./", length=1.0, output = None, skip_zeros=True, show_plot=False):
		"""
		Reads a file containing several timesteps with average velocity profile across two parallel plates.
		Function will take average of all timesteps in order to get good statistics.
		"""
		temperatures = pylab.loadtxt(state_path+"/statistics/linear_temperature.txt")
		num_bins = len(temperatures)   # Number of elements on first line equals number of bins
		print "Plotting linear temperature profile with "+str(num_bins)+" bins."

		x_axis = linspace(0,length,num_bins) # Create x_axis with possible real channel height

		if skip_zeros:
			self.plot(x=x_axis[temperatures>0], y=temperatures[temperatures>0], x_axis_label="x", y_axis_label="Temperature [K]", filename=output)
		else:
			self.plot(x=x_axis, y=temperatures, x_axis_label="x", y_axis_label="Temperature [K]", filename=output)
		
		if show_plot: plt.show()

	def plot_linear_pressure_profile(self, state_path="./", length=1.0, output = None, skip_zeros=True, show_plot=False):
		"""
		Reads a file containing several timesteps with average velocity profile across two parallel plates.
		Function will take average of all timesteps in order to get good statistics.
		"""
		pressures = pylab.loadtxt(state_path+"/statistics/linear_pressure.txt")
		num_bins = len(pressures)   # Number of elements on first line equals number of bins
		print "Plotting linear pressure profile with "+str(num_bins)+" bins."

		x_axis = linspace(0,length,num_bins) # Create x_axis with possible real channel height

		if skip_zeros:
			self.plot(x=x_axis[pressures>0], y=pressures[pressures>0], x_axis_label="x", y_axis_label="Pressure [Pa]", filename=output)
		else:
			self.plot(x=x_axis, y=pressures, x_axis_label="x", y_axis_label="Pressure [Pa]", filename=output)
			
		if show_plot: plt.show()

	def plot_permeability(self, state_path="./", output = None, show_plot=False):
		"""
		Reads a file containing several timesteps with average velocity profile across a cylinder.
		Function will take average of all timesteps in order to get good statistics.
		"""

		data = pylab.loadtxt(state_path+"/statistics/permeability.txt")
		num_timesteps = len(data) # Number of lines equals number of timesteps
		time = data[:,0]   # Number of elements on first line equals number of bins
		permeabilities = data[:,1]   # Number of elements on first line equals number of bins
		
		print "Plotting time vs permeability during "+str(num_timesteps)+" timesteps."

		self.plot(x=time, y=permeabilities, x_axis_label="t", y_axis_label="k", filename=output)
		if show_plot: plt.show()