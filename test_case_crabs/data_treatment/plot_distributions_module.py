# --------------------------------------------------------------------------------------------------------------
# BEAM DISTRIBUTION PLOT
# ANDREA SANTAMARIA
# LAST MODFIIED: 26/08/2014
# This script plots the generated distribution
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# Import libraries and introduce abbreviations
# --------------------------------------------------------------------------------------------------------------
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import ticker
from math import *
from collections import Counter
import pylab as P
import scipy.stats
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
from matplotlib.ticker import LogLocator, LogFormatter 
from matplotlib.colors import LogNorm
from matplotlib.colors import BoundaryNorm
from matplotlib.colorbar import ColorbarBase
plt.ion()



# --------------------------------------------------------------------------------------------------------------
# Text characteristics of the plot
# --------------------------------------------------------------------------------------------------------------
params 		= {	'backend': 'pdf',
	  		'font.size':20,	
          	'axes.labelsize': 20,
          	'legend.fontsize': 16,
          	'xtick.labelsize': 20,
          	'ytick.labelsize': 20,
	  		'text.usetex':True,
        	 } 
# --------------------------------------------------------------------------------------------------------------
# Activate LaTeX 
# --------------------------------------------------------------------------------------------------------------        
rc('text.latex', preamble=r'\usepackage{cmbright}')
matplotlib.rcParams.update(params)

def beam_histogram(n_part, variable, string, title):

	fig 		= figure()
	ax 			= fig.add_subplot(1, 1, 1)

	mu 			= np.mean(variable)
	sigma 		= np.std(variable)

	# To put it in units of sigma 
	# --------------------------------------------------------------------------------------------------------------
	var_sigma	= variable/sigma
	var_max 	= max(var_sigma)
	var_min 	= min(var_sigma)

	# --------------------------------------------------------------------------------------------------------------
	# Compute the best number of bins 
	# --------------------------------------------------------------------------------------------------------------

	# Interquartile range (IQR)
	# --------------------------------------------------------------------------------------------------------------
	IQR 		= np.percentile(var_sigma, 0.75) - np.percentile(var_sigma, 0.25)

	# Bin size following the Freedman Diaconis rule
	# --------------------------------------------------------------------------------------------------------------
	bin_size 	= 2 * IQR * n_part**(-1.0/3)
	print 'Bin size = %s'% bin_size

	# Number of bins
	# --------------------------------------------------------------------------------------------------------------
	nbins 		= int(var_max - var_min/ (bin_size))
	print 'Number of bins= %s'% nbins

	# --------------------------------------------------------------------------------------------------------------
	# Plot Histogram of Samples
	# --------------------------------------------------------------------------------------------------------------
	n, bins, patches = P.hist((var_sigma), bins = nbins-0.8*nbins, normed = True, histtype = 'stepfilled', label = 'Histogram of the samples')
	P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

	# --------------------------------------------------------------------------------------------------------------
	# Plot the Probability Density Function
	# --------------------------------------------------------------------------------------------------------------
	mu_new			= np.mean(var_sigma)
	sigma_new		= np.std(var_sigma)
	y 				= P.normpdf(bins, mu_new, sigma_new)
	P.plot(bins, y, 'k--', linewidth=1.5, label = 'Probability Density Function')

	density = gaussian_kde(var_sigma)  # your data
	xgrid = np.linspace(var_min, var_max, n_part)   
	plt.plot(xgrid, density(xgrid),'k--', color = 'red', linewidth=1.5, label = 'Kernel Density Estimation Method')

	# Horizontal axis
	# --------------------------------------------------------------------------------------------------------------
	plt.xlabel(r'%s [$\sigma$]'%string)
	plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0))
	plt.xlim(var_min, var_max)
	# plt.grid(b = True, which = 'both', axis = 'both', linestyle = '--')

	# Vertical axis
	# --------------------------------------------------------------------------------------------------------------
	plt.ylabel(r'Normalized beam profile')
	# ax.set_yscale('log')
	# plt.ylim(0, max(pdf_fitted))
	# plt.grid(b = True, which = 'minor', axis = 'both', linestyle = '--')
	
	# Title, legends and annotations
	# --------------------------------------------------------------------------------------------------------------
	# title 		= 'Distribution in : particles = %.0f, $\mu$ = %E, $\sigma$ = %E'%(n_part,mu,sigma)
	# plt.title(title)
	# plt.legend(loc = 'upper left')
	# plt.text((var_max/2) + (var_max/14), max(variable)/2 + max(variable)/3.0, r'Particles = %.0f'%n_part)
	# plt.text((var_max/2) + (var_max/14), max(variable)/2 + max(variable)/3.7 , r'$\mu$ = %E'%mu)
	# plt.text((var_max/2) + (var_max/14), max(variable)/2 + max(variable)/4.5, r'$\sigma$ = %E'%sigma)   

	# Save to a File
	filename = '%s'%title
	plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)	


	

def beam_profile(n_part, variable, string, title):

	fig 		= figure()
	ax 			= fig.add_subplot(1, 1, 1)

	mu 			= np.mean(variable)
	sigma 		= np.std(variable)

	# To put it in units of sigma 
	# --------------------------------------------------------------------------------------------------------------
	var_sigma	= variable/sigma
	var_max 	= max(var_sigma)
	var_min 	= min(var_sigma)

	mu_new 		= np.mean(var_sigma)
	sigma_new 	= np.std(var_sigma)

	# --------------------------------------------------------------------------------------------------------------
	# Compute the best number of bins 
	# --------------------------------------------------------------------------------------------------------------

	# Interquartile range (IQR)
	# --------------------------------------------------------------------------------------------------------------
	IQR 		= np.percentile(var_sigma, 0.75) - np.percentile(var_sigma, 0.25)

	# Bin size following the Freedman Diaconis rule
	# --------------------------------------------------------------------------------------------------------------
	bin_size 	= 2 * IQR * n_part**(-1.0/3)
	print 'Bin size = %s'% bin_size

	# Number of bins
	# --------------------------------------------------------------------------------------------------------------
	nbins 		= int(var_max - var_min/ (bin_size))
	print 'Number of bins= %s'% nbins

	# --------------------------------------------------------------------------------------------------------------
	# Plot Histogram of Samples
	# --------------------------------------------------------------------------------------------------------------
	n, bins, patches = P.hist((var_sigma), bins = nbins-0.0*nbins, normed = False, histtype = 'step', label = 'Beam Profile')
	P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)


	# Horizontal axis
	# --------------------------------------------------------------------------------------------------------------
	plt.xlabel(r'%s [$\sigma$]'%string)
	plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0))
	plt.xlim(var_min, var_max)
	# plt.grid(b = True, which = 'both', axis = 'both', linestyle = '--')

	# Vertical axis
	# --------------------------------------------------------------------------------------------------------------
	plt.ylabel(r'Number of Particles')
	# ax.set_yscale('log')
	# plt.ylim(0, max(pdf_fitted))
	# plt.grid(b = True, which = 'minor', axis = 'both', linestyle = '--')
	
	# Title, legends and annotations
	# --------------------------------------------------------------------------------------------------------------
	# title 		= 'Distribution in : particles = %.0f, $\mu$ = %E, $\sigma$ = %E'%(n_part,mu,sigma)
	# plt.title(title)
	# plt.legend(loc = 'upper left')
	# plt.text((var_max/2) + (var_max/14), (n_part/20)/(sqrt(2*pi)*sigma_new), r'Particles = %.0f'%n_part)
	# plt.text((var_max/2) + (var_max/14), (n_part/22)/(sqrt(2*pi)*sigma_new), r'$\sigma$ = %E'%sigma)


	# # Plot exactly what's in the file (no normalisation)
	# # --------------------------------------------------------------------------------------------------------------
	# c 			= Counter(var_sigma)
	# x_values 	= c.keys()
	# y_values	= c.values()
	# ax.bar(x_values, y_values, log = False, width = bin_size, align = 'center')
	# plt.xlabel(r'%s [$\sigma$]'%string)
	# plt.ylabel(r'Number of Particles')
	# plt.text((var_max/2) + (var_max/14), max(variable_2)/2 + max(variable_2)/3.1, r'Particles = %.0f'%n_part)
	# 
	
	# Save to a File
	filename = '%s'%title
	plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)	



def beam_density(variable_1, variable_2, string_1, string_2, n_part, title):

	# Define a function to make the ellipses
	def ellipse(ra,rb,ang,x0,y0,Nb=100):
	    xpos,ypos=x0,y0
	    radm,radn=ra,rb
	    an=ang
	    co,si=np.cos(an),np.sin(an)
	    the=linspace(0,2*np.pi,Nb)
	    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
	    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
	    return X,Y
	# --------------------------------------------------------------------------------------------------------------
	# VARIABLE 1
	# --------------------------------------------------------------------------------------------------------------
	mu_1 		= np.mean(variable_1)
	sigma_1 	= np.std(variable_1)

	var_sigma_1	= variable_1 / sigma_1
	var_max_1 	= max(var_sigma_1)
	var_min_1 	= min(var_sigma_1)

	mu_new_1 	= np.mean(var_sigma_1)
	sigma_new_1 = np.std(var_sigma_1)

	# Compute the best number of bins 
	# --------------------------------------------------------------------------------------------------------------

	# Interquartile range (IQR)
	# --------------------------------------------------------------------------------------------------------------
	IQR 		= np.percentile(var_sigma_1, 0.75) - np.percentile(var_sigma_1, 0.25)

	# Bin size following the Freedman Diaconis rule
	# --------------------------------------------------------------------------------------------------------------
	bin_size_1 	= 2 * IQR * n_part**(-1.0/3)
	print 'Bin size 1 = %s'% bin_size_1

	# Number of bins
	# --------------------------------------------------------------------------------------------------------------
	#Set up the histogram bins
	nbins_1 	= int((var_max_1 - var_min_1) / bin_size_1)
	print 'Number of bins= %s'% nbins_1


	# --------------------------------------------------------------------------------------------------------------
	# VARIABLE 2
	# --------------------------------------------------------------------------------------------------------------
	mu_2 		= np.mean(variable_2)
	sigma_2 	= np.std(variable_2)

	var_sigma_2	= variable_2 / sigma_2
	var_max_2 	= max(var_sigma_2)
	var_min_2 	= min(var_sigma_2)

	mu_new_2 	= np.mean(var_sigma_2)
	sigma_new_2 = np.std(var_sigma_2)

	# Compute the best number of bins 
	# --------------------------------------------------------------------------------------------------------------
	
	# Interquartile range (IQR)
	# --------------------------------------------------------------------------------------------------------------
	IQR 		= np.percentile(var_sigma_2, 0.75) - np.percentile(var_sigma_2, 0.25)

	# Bin size following the Freedman Diaconis rule
	# --------------------------------------------------------------------------------------------------------------
	bin_size_2 	= 2 * IQR * n_part**(-1.0/3)
	print 'Bin size 1 = %s'% bin_size_2

	# Number of bins
	# --------------------------------------------------------------------------------------------------------------
	nbins_2 	= int((var_max_2 - var_min_2) / bin_size_2)
	print 'Number of bins= %s'% nbins_2

	# Set up default x and y limits
	xlims = [var_min_1, var_max_1]
	ylims = [var_min_2, var_max_2]

	# Set up your x and y labels
	xlabel = r'%s'%string_1
	ylabel = r'%s'%string_2

	# Define the locations for the axes
	left, width = 0.12, 0.55
	bottom, height = 0.12, 0.55
	bottom_h = left_h = left+width+0.02
	 
	# Set up the geometry of the three plots
	rect_temperature = [left, bottom, width, height] # dimensions of temp plot
	rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
	rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
	 
	# Set up the size of the figure
	fig = plt.figure(4, figsize=(9.5,9))
		 
	# Make the three plots
	axTemperature 		= plt.axes(rect_temperature) # temperature plot
	axHistx 			= plt.axes(rect_histx) # x histogram
	axHisty 			= plt.axes(rect_histy) # y histogram
	 
	# Remove the inner axes numbers of the histograms
	nullfmt 			= NullFormatter()
	axHistx.xaxis.set_major_formatter(nullfmt)
	axHisty.yaxis.set_major_formatter(nullfmt)
	 
	# Find the min/max of the data

	 
	# Make the 'main' temperature plot
	# Define the number of bins
	 
	xbins 				= linspace(start = var_min_1, stop = var_max_1, num = nbins_1)
	ybins 				= linspace(start = var_min_2, stop = var_max_2, num = nbins_2)
	xcenter 			= (xbins[0:-1]+xbins[1:])/2.0
	ycenter 			= (ybins[0:-1]+ybins[1:])/2.0
	aspectratio 		= 1.0*(var_max_1 - 0)/(1.0*var_max_2 - 0)
	
	nx=int(nbins_1 - 0*(nbins_1))
	ny=int(nbins_2 - 0*(nbins_2))

	 
	H, xedges,yedges 	= np.histogram2d(var_sigma_1, var_sigma_2, bins=[nx,ny])
	X 					= xcenter
	Y 					= ycenter
	Z 					= H
	 
	# Plot the temperature data
	cax 				= (axTemperature.imshow(H, extent=[var_min_1,var_max_1,var_min_2,var_max_2],
	       					interpolation='nearest', origin='lower',aspect=aspectratio))
	 
	# Plot the temperature plot contours
	contourcolor 		= 'white'
	xcenter 			= mu_new_1
	ycenter 			= mu_new_2
	ra 					= sigma_new_1
	rb 					= sigma_new_2
	ang 				= 0
	 
	X,Y 				= ellipse(ra,rb,ang,xcenter,ycenter)
	axTemperature.plot(X,Y,"k:",ms=1,linewidth=2.0)
	axTemperature.annotate('$1\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
	                       textcoords='offset points', horizontalalignment='right',
	                       verticalalignment='bottom',fontsize=25)
	 
	X,Y=ellipse(2*ra,2*rb,ang,xcenter,ycenter)
	axTemperature.plot(X,Y,"k:",color = contourcolor,ms=1,linewidth=2.0)
	axTemperature.annotate('$2\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
	                       textcoords='offset points',horizontalalignment='right',
	                       verticalalignment='bottom',fontsize=25, color = contourcolor)
	 
	X,Y=ellipse(3*ra,3*rb,ang,xcenter,ycenter)
	axTemperature.plot(X,Y,"k:",color = contourcolor, ms=1,linewidth=2.0)
	axTemperature.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
	                       textcoords='offset points',horizontalalignment='right',
	                       verticalalignment='bottom',fontsize=25, color = contourcolor)
	 
	#Plot the axes labels
	axTemperature.set_xlabel(xlabel,fontsize=25)
	axTemperature.set_ylabel(ylabel,fontsize=25)
	 
	#Make the tickmarks pretty
	ticklabels = axTemperature.get_xticklabels()
	for label in ticklabels:
	    label.set_fontsize(18)
	    label.set_family('serif')
	 
	ticklabels = axTemperature.get_yticklabels()
	for label in ticklabels:
	    label.set_fontsize(18)
	    label.set_family('serif')
	 
	#Set up the plot limits
	axTemperature.set_xlim(xlims)
	axTemperature.set_ylim(ylims)
	 
	#Plot the histograms
	axHistx.hist(var_sigma_1, bins = nbins_1-0.0*nbins_1, color = 'blue', normed = True, histtype = 'stepfilled')
	axHisty.hist(var_sigma_2, bins = nbins_2-0.0*nbins_2, orientation = 'horizontal', color = 'red', normed = True, histtype = 'stepfilled')
	 
	#Set up the histogram limits
	# axHistx.set_xlim(xlims)
	# axHisty.set_ylim(ylims)
	 
	#Make the tickmarks pretty
	ticklabels = axHistx.get_yticklabels()
	for label in ticklabels:
	    label.set_fontsize(12)
	    label.set_family('serif')
	 
	#Make the tickmarks pretty
	ticklabels = axHisty.get_xticklabels()
	for label in ticklabels:
	    label.set_fontsize(12)
	    label.set_family('serif')
	 
	#Cool trick that changes the number of tickmarks for the histogram axes
	axHisty.xaxis.set_major_locator(MaxNLocator(4))
	axHistx.yaxis.set_major_locator(MaxNLocator(4))

	#Show the plot
	# plt.draw()
	 
	# Save to a File
	filename = '%s'%title
	plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)	


def beam_2dhist_comparison(variable_1, variable_2, variable_3, variable_4, string_1, string_2, n_part, title):

	# fig 	= figure()
	fig 	= plt.figure(2, figsize=(16,14))
	ax1 	= fig.add_subplot(121)
	ax2 	= fig.add_subplot(122)

	# --------------------------------------------------------------------------------------------------------------
	# VARIABLE 1
	# --------------------------------------------------------------------------------------------------------------
	mu_1 		= np.mean(variable_1)
	sigma_1 	= np.std(variable_1)

	var_max_1	= max(variable_1)
	var_min_1	= min(variable_1)

	# Compute the best number of bins 
	# --------------------------------------------------------------------------------------------------------------

	# Interquartile range (IQR)
	# --------------------------------------------------------------------------------------------------------------
	IQR 		= np.percentile(variable_1, 0.75) - np.percentile(variable_1, 0.25)

	# Bin size following the Freedman Diaconis rule
	# --------------------------------------------------------------------------------------------------------------
	bin_size_1 	= 2 * IQR * n_part**(-1.0/3)
	print 'Bin size 1 = %s'% bin_size_1

	# Number of bins
	# --------------------------------------------------------------------------------------------------------------
	#Set up the histogram bins
	nbins_1 	= int((var_max_1 - var_min_1) / bin_size_1)
	print 'Number of bins= %s'% nbins_1


	# --------------------------------------------------------------------------------------------------------------
	# VARIABLE 2
	# --------------------------------------------------------------------------------------------------------------
	mu_2 		= np.mean(variable_2)
	sigma_2 	= np.std(variable_2)

	var_max_2 	= max(variable_2)
	var_min_2 	= min(variable_2)

	# Compute the best number of bins 
	# --------------------------------------------------------------------------------------------------------------
	
	# Interquartile range (IQR)
	# --------------------------------------------------------------------------------------------------------------
	IQR 		= np.percentile(variable_2, 0.75) - np.percentile(variable_2, 0.25)

	# Bin size following the Freedman Diaconis rule
	# --------------------------------------------------------------------------------------------------------------
	bin_size_2 	= 2 * IQR * n_part**(-1.0/3)
	print 'Bin size 2 = %s'% bin_size_2

	# Number of bins
	# --------------------------------------------------------------------------------------------------------------
	nbins_2 	= int((var_max_2 - var_min_2) / bin_size_2)
	print 'Number of bins= %s'% nbins_2


	# Set up default x and y limits
	# --------------------------------------------------------------------------------------------------------------
	xlims_1 = [var_min_1, var_max_1]
	ylims_1 = [var_min_2, var_max_2]


	# --------------------------------------------------------------------------------------------------------------
	# VARIABLE 3
	# --------------------------------------------------------------------------------------------------------------
	mu_3 		= np.mean(variable_3)
	sigma_3 	= np.std(variable_3)

	var_max_3	= max(variable_3)
	var_min_3	= min(variable_3)

	# Compute the best number of bins 
	# --------------------------------------------------------------------------------------------------------------

	# Interquartile range (IQR)
	# --------------------------------------------------------------------------------------------------------------
	IQR 		= np.percentile(variable_3, 0.75) - np.percentile(variable_3, 0.25)

	# Bin size following the Freedman Diaconis rule
	# --------------------------------------------------------------------------------------------------------------
	bin_size_3 	= 2 * IQR * n_part**(-1.0/3)
	print 'Bin size 3 = %s'% bin_size_3

	# Number of bins
	# --------------------------------------------------------------------------------------------------------------
	#Set up the histogram bins
	nbins_3 	= int((var_max_3 - var_min_3) / bin_size_3)
	print 'Number of bins= %s'% nbins_3


	# --------------------------------------------------------------------------------------------------------------
	# VARIABLE 4
	# --------------------------------------------------------------------------------------------------------------
	mu_4 		= np.mean(variable_4)
	sigma_4 	= np.std(variable_4)

	var_max_4	= max(variable_4)
	var_min_4	= min(variable_4)

	# Compute the best number of bins 
	# --------------------------------------------------------------------------------------------------------------

	# Interquartile range (IQR)
	# --------------------------------------------------------------------------------------------------------------
	IQR 		= np.percentile(variable_4, 0.75) - np.percentile(variable_4, 0.25)

	# Bin size following the Freedman Diaconis rule
	# --------------------------------------------------------------------------------------------------------------
	bin_size_4 	= 2 * IQR * n_part**(-1.0/3)
	print 'Bin size 4 = %s'% bin_size_4

	# Number of bins
	# --------------------------------------------------------------------------------------------------------------
	#Set up the histogram bins
	nbins_4 	= int((var_max_4 - var_min_4) / bin_size_4)
	print 'Number of bins= %s'% nbins_4

	# Set up default x and y limits
	# --------------------------------------------------------------------------------------------------------------
	xlims_2 = [var_min_3, var_max_3]
	ylims_2 = [var_min_4, var_max_4]



	#FIRST PLOT (INITIAL CONDITIONS)
	# --------------------------------------------------------------------------------------------------------------
	
	xbins_1 			= linspace(start = var_min_1, stop = var_max_1, num = sqrt(nbins_1) + 0.08*nbins_1)
	ybins_1 			= linspace(start = var_min_2, stop = var_max_2, num = sqrt(nbins_2) + 0.08*nbins_2)

	xcenter 			= (xbins_1[0:-1]+xbins_1[1:])/2.0
	ycenter 			= (ybins_1[0:-1]+ybins_1[1:])/2.0
	aspectratio 		= 1.0*(var_max_1 - 0)/(1.0*var_max_2 - 0)

	H, xedges,yedges 	= np.histogram2d(variable_2, variable_1, bins=(ybins_1, xbins_1))
	X 					= xcenter
	Y 					= ycenter
	Z 					= H
	cax 				= (ax1.imshow(H, extent=[var_min_1,var_max_1,var_min_2,var_max_2],
       					interpolation='nearest', origin='lower',aspect=aspectratio))
	ax1.set_xlabel(r'%s'%string_1)
	ax1.set_ylabel(r'%s'%string_2)

	#Set up the plot limits
	ax1.set_xlim(xlims_1)
	ax1.set_ylim(ylims_1)

	#cbar = plt.colorbar(cax, ax= ax1, shrink=.4)

	ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
	
	#SECOND PLOT (FINAL CONDITIONS)
	# --------------------------------------------------------------------------------------------------------------
	xbins_2 			= linspace(start = var_min_3, stop = var_max_3, num = sqrt(nbins_3) + 0.08*nbins_3)
	ybins_2 			= linspace(start = var_min_4, stop = var_max_4, num = sqrt(nbins_4) + 0.08*nbins_4)

	xcenter 			= (xbins_2[0:-1]+xbins_2[1:])/2.0
	ycenter 			= (ybins_2[0:-1]+ybins_2[1:])/2.0
	aspectratio 		= 1.0*(var_max_3 - 0)/(1.0*var_max_4 - 0)

	cmap 	= plt.cm.jet
	bounds 	= np.linspace(0,0.9,1)
	norm 	= BoundaryNorm(bounds, cmap.N)

	H, xedges,yedges 	= np.histogram2d(variable_4, variable_3, bins=(ybins_2, xbins_2))
	cax 				= (ax2.imshow(H, cmap = cmap, norm = norm  ,extent=[var_min_3,var_max_3,var_min_4,var_max_4],
       					interpolation='nearest', origin='lower',aspect=aspectratio))

	ax2.set_xlabel(r'%s'%string_1)
	ax2.set_ylabel(r'%s'%string_2)

	#Set up the plot limits
	ax2.set_xlim(xlims_2)
	ax2.set_ylim(ylims_2)
	ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

	#ax3 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
	#cbar = plt.colorbar(ax3, shrink=.4, cmap = cmap, norm = norm, ticks=bounds, boundaries=bounds)
	# ax3.set_label('Density (arbitrary units)')


	#Show the plot
	# plt.draw()
	 
	# Save to a File
	filename = '%s'%title
	plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)	


def beam_scatter(variable_1, variable_2, string_1, string_2):

	fig 		= figure()
	ax 			= fig.add_subplot(1, 1, 1)


	a = np.vstack([variable_1,variable_2])
	t = gaussian_kde(a)(a)						# This uses multiple normal distributions (gaussian functions) to estimate 
												# the probability density function. 
												# kde = kernel density estimation: is a non-parametric way to estimate the 
												# probability density function of a random variable.

	ax.scatter(variable_1, variable_2, c = t)
	ax.set_xlabel(r'%s'%string_1)
	ax.set_ylabel(r'%s'%string_2)
	ax.set_xlim([min(variable_1) + (min(variable_1)), max(variable_1) + (max(variable_1))])
	ax.set_ylim([min(variable_2) + (min(variable_2)), max(variable_2) + (max(variable_2))])
	ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
	ax.grid(b=True, which='major',linestyle='--')

	

def beam_scatter_three(variable_1, variable_2, variable_3, string_1, string_2, string_3):

	fig 		= figure()
	ax 			= fig.add_subplot(111, projection='3d')

	ax.scatter(variable_1, variable_2, variable_3)
	ax.set_xlabel('%s'%string_1)
	ax.set_ylabel('%s'%string_2)
	ax.set_zlabel('%s'%string_3)
	ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

	

def beam_scatter_comparison(variable_1, variable_2, variable_3, variable_4, string_1, string_2, title):

	fig 	= figure()
	ax1 	= fig.add_subplot(121)
	ax2 	= fig.add_subplot(122)

	#FIRST PLOT (INITIAL CONDITIONS)
	# --------------------------------------------------------------------------------------------------------------

	a 		= np.vstack([variable_1,variable_2])
	t		= gaussian_kde(a)(a)

	ax1.scatter(variable_1, variable_2, c = t)
	ax1.set_xlabel(r'%s'%string_1)
	ax1.set_ylabel(r'%s'%string_2)
	ax1.set_xlim([min(variable_1) + (min(variable_1)), max(variable_1) + (max(variable_1))])
	ax1.set_ylim([min(variable_2) + (min(variable_2)), max(variable_2) + (max(variable_2))])
	ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
	ax1.grid(b=True, which='major',linestyle='--')
	
	#SECOND PLOT (FINAL CONDITIONS)
	# --------------------------------------------------------------------------------------------------------------

	b = np.vstack([variable_3, variable_4])
	m = gaussian_kde(b)(b)

	ax2.scatter(variable_3, variable_4, c = m)
	ax2.set_xlabel(r'%s'%string_1)
	ax2.set_ylabel(r'%s'%string_2)
	ax2.set_xlim([min(variable_3) + (min(variable_3)), max(variable_3) + (max(variable_3))])
	ax2.set_ylim([min(variable_4) + (min(variable_4)), max(variable_4) + (max(variable_4))])
	ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
	ax2.grid(b=True, which='major',linestyle='--')

	#COMMON TITLE
	# --------------------------------------------------------------------------------------------------------------
	pyplot.suptitle('Initial and final particle distributions', fontsize=30)

	# Save to a File
	filename = '%s'%title
	plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)	



def beam_scatter_comparison_three(variable_1, variable_2, variable_3, variable_4, variable_5, variable_6 ):

	fig = plt.figure()
	ax1 = fig.add_subplot(121, projection='3d')
	ax2 = fig.add_subplot(122, projection='3d')

	ax1.scatter(variable_3, variable_1, variable_2)
	ax1.set_xlabel(r'z (m)')
	ax1.set_ylabel(r'x (m)')
	ax1.set_zlabel(r'y (m)')
	ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))


	ax2.scatter(variable_6, variable_4, variable_5)
	ax2.set_xlabel(r'z (m)')
	ax2.set_ylabel(r'x (m)')
	ax2.set_zlabel(r'y (m)')
	ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

	#COMMON TITLE
	# --------------------------------------------------------------------------------------------------------------
	pyplot.suptitle('Initial and final particle distributions', fontsize=30)






