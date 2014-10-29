# --------------------------------------------------------------------------------------------------------------
# BEAM DISTRIBUTION COMPARISON AFTER SIMULATION
# ANDREA SANTAMARIA
# LAST MODFIIED: 02/09/2014
# 
# This script allows you to:
# 	- Treat several output files from SixTrack and compare them
# --------------------------------------------------------------------------------------------------------------

import 	sys
import 	os
import 	matplotlib
from 	matplotlib import pyplot as plt

#sys.path.append("/afs/cern.ch/work/a/ansantam/private/crab_cavity_simulations/python_modules/")

from 	data_treatment_module import *
from 	plot_distributions_module import *
from 	plot_absorptions_module import *


# --------------------------------------------------------------------------------------------------------------
# Load each parameter from the initial distribution to use it in the plots
# --------------------------------------------------------------------------------------------------------------
(x, xp, y, yp, z, e) 		= np.loadtxt(r'../output/dist0.dat', unpack=True)
(xn, xpn, yn, ypn, zn, en) 	= np.loadtxt(r'../output/distn.dat', unpack=True)

 # --------------------------------------------------------------------------------------------------------------
 # Check that the number of lines = number of particles simulated and lost
 # --------------------------------------------------------------------------------------------------------------
n_part 	= len(x)
print 'Number of particles initial distribution 	= %s'% n_part

n_part_final 	= len(xn)
print 'Number of particles final distribution 		= %s'% n_part_final

lost 		= n_part - n_part_final 
print 'Number of particles lost 			= %s'% lost

percentage_lost = (lost * 100) / n_part
print 'Percentage of particles lost 			= %s'% percentage_lost


# --------------------------------------------------------------------------------------------------------------
# Plot the initial and final distributions for comparison
# --------------------------------------------------------------------------------------------------------------

# 3D overview of the bunches
# --------------------------------------------------------------------------------------------------------------
#beam_scatter_comparison_three(x, y, z, xn, yn, zn)

# X coordinate
# --------------------------------------------------------------------------------------------------------------
beam_scatter_comparison(x, xp, xn, xpn, 'x', 'xp','x_scatter')

# Y coordinate
# --------------------------------------------------------------------------------------------------------------
beam_scatter_comparison(y, yp, yn, ypn, 'y', 'yp', 'y_scatter')

# Z coordinate
# --------------------------------------------------------------------------------------------------------------
#beam_scatter_comparison(z, e, zn, en, 'z', 'e', 'ze_scatter')
beam_scatter_comparison(z, x, zn, xn, 'z', 'x', 'zx_scatter')
beam_scatter_comparison(z, y, zn, yn, 'z', 'y', 'zy_scatter')

# Histograms
# --------------------------------------------------------------------------------------------------------------
#beam_histogram(n_part, xn, 'x', 'x_histogram_final')
#beam_profile(n_part, yn, 'y', 'x_profile_final')

#beam_histogram(n_part, x, 'x', 'y_histogram')
#beam_profile(n_part, y, 'y', 'y_profile')
#beam_density(x, xp, 'x', 'xp', n_part,'x_density')
#beam_density(y, yp, 'y', 'yp', n_part,'y_density')

#beam_2dhist_comparison(x, xp, xn, xpn, 'x', "$x'$", n_part, 'x_comparison')

# Absorptions
#--------------------------------------------------------------------------------------------------------------
#convert_to_csv("../output/all_absorptions.dat", "abs.csv")
#convert_to_csv("../output/LPI_test.s", "ap.csv")

#plot_absorptions("abs.csv", "ap.csv")

plt.show()
