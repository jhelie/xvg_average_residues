#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'xvg_average_residues_simple', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_average_op
**********************************************

[ DESCRIPTION ]
 
This script calculate the average of clustering data contained in several xvg files.

It calculates the avg and (unbiasd) std dev and can deal with NaN.

NB:
the script may give out a warning 'return np.mean(x,axis)/factor', it's ok. it's just
scipy warning us that there were only nans on a row, the result will be a nan as we
expect (see this thread: https://github.com/scipy/scipy/issues/2898).

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file(s)
-o	residues_avg	: name of outptut file
--membrane		: 'AM_zCter','AM_zNter','SMa' or 'SMz'
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs='+', dest='xvgfilenames', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["residues_avg"], help=argparse.SUPPRESS)
parser.add_argument('--membrane', dest='membrane', choices=['AM_zCter','AM_zNter','SMa','SMz'], default='not specified', help=argparse.SUPPRESS, required=True)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.output_file = args.output_file[0]
args.comments = args.comments[0].split(',')

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

if len(args.xvgfilenames) == 1:
	print "Error: only 1 data file specified."
	sys.exit(1)
	
for f in args.xvgfilenames:
	if not os.path.isfile(f):
		print "Error: file " + str(f) + " not found."
		sys.exit(1)

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
	
	#NB
	#Asymmetric membranes:
	#  - upper = zwitterionic leaflet
	#  - lower = anionic leaflet
	#
	#Symmetric membranes:
	#  - upper = Nter leaflet
	#  - lower = Cter leaflet
	
	global nb_rows
	global nb_cols
	global weights
	global distance
	global data_residues_basic
	global data_residues_upper_POPE
	global data_clustering_upper_POPS
	global data_clustering_lower_POPC
	global data_clustering_lower_POPE
	global data_clustering_lower_POPS

	nb_rows = 0
	nb_cols = 0
	weights = np.ones(len(args.xvgfilenames))
		
	for f_index in range(0,len(args.xvgfilenames)):
		#display progress
		progress = '\r -reading file ' + str(f_index+1) + '/' + str(len(args.xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		
		#get file content
		filename = args.xvgfilenames[f_index]
		with open(filename) as f:
			lines = f.readlines()
		
		#determine legends and nb of lines to skip
		tmp_nb_rows_to_skip = 0
		for l_index in range(0,len(lines)):
			line = lines[l_index]
			if line[0] in args.comments:
				tmp_nb_rows_to_skip += 1
				if "weight" in line:
					if "-> weight = " in line:
						weights[f_index] = float(line.split("-> weight = ")[1])
						if weights[f_index] < 0:
							print "\nError: the weight in file " + str(filename) + " should be a positive number."
							print " -> " + str(line)
							sys.exit(1)
					else:
						print "\nWarning: keyword 'weight' found in the comments of file " + str(filename) + ", but weight not read in as the format '-> weight = ' wasn't found."
		
		#get data
		tmp_data = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
		
		#check that each file has the same number of data rows
		if f_index == 0:
			nb_rows = np.shape(tmp_data)[0]
			distance = np.zeros((nb_rows, 1))													#distance from cluster
			data_clustering_upper_POPC = np.zeros((nb_rows, len(args.xvgfilenames)))			#upper POPC %
			data_clustering_upper_POPE = np.zeros((nb_rows, len(args.xvgfilenames)))			#upper POPE %
			data_clustering_upper_POPS = np.zeros((nb_rows, len(args.xvgfilenames)))			#upper POPS %
			data_clustering_lower_POPC = np.zeros((nb_rows, len(args.xvgfilenames)))			#lower POPC %
			data_clustering_lower_POPE = np.zeros((nb_rows, len(args.xvgfilenames)))			#lower POPE %
			data_clustering_lower_POPS = np.zeros((nb_rows, len(args.xvgfilenames)))			#lower POPS %
		else:
			if np.shape(tmp_data)[0] != nb_rows:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[0]) + " data rows, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_rows) + " data rows."
				sys.exit(1)
		#check that each file has the same number of columns
		if f_index == 0:
			nb_cols = np.shape(tmp_data)[1]
		else:
			if np.shape(tmp_data)[1] != nb_cols:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[1]) + " data columns, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_cols) + " data columns."
				sys.exit(1)
		#check that each file has the same first column
		if f_index == 0:
			distance[:,0] = tmp_data[:,0]
		else:
			if not np.array_equal(tmp_data[:,0],distance[:,0]):
				print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
				sys.exit(1)
		
		#store data
		if args.membrane == "AM_zCter":
			data_clustering_upper_POPC[:,f_index] = tmp_data[:,3]
			data_clustering_upper_POPE[:,f_index] = tmp_data[:,4]
			data_clustering_lower_POPC[:,f_index] = tmp_data[:,8]
			data_clustering_lower_POPE[:,f_index] = tmp_data[:,9]
			data_clustering_lower_POPS[:,f_index] = tmp_data[:,10]
		elif args.membrane == "AM_zNter":
			data_clustering_upper_POPC[:,f_index] = tmp_data[:,9]
			data_clustering_upper_POPE[:,f_index] = tmp_data[:,10]
			data_clustering_lower_POPC[:,f_index] = tmp_data[:,4]
			data_clustering_lower_POPE[:,f_index] = tmp_data[:,5]
			data_clustering_lower_POPS[:,f_index] = tmp_data[:,6]
		elif args.membrane == "SMa":
			data_clustering_upper_POPC[:,f_index] = tmp_data[:,10]
			data_clustering_upper_POPE[:,f_index] = tmp_data[:,11]
			data_clustering_upper_POPS[:,f_index] = tmp_data[:,12]
			data_clustering_lower_POPC[:,f_index] = tmp_data[:,4]
			data_clustering_lower_POPE[:,f_index] = tmp_data[:,5]
			data_clustering_lower_POPS[:,f_index] = tmp_data[:,6]
		elif args.membrane == "SMz":
			data_clustering_upper_POPC[:,f_index] = tmp_data[:,7]
			data_clustering_upper_POPE[:,f_index] = tmp_data[:,8]
			data_clustering_lower_POPC[:,f_index] = tmp_data[:,3]
			data_clustering_lower_POPE[:,f_index] = tmp_data[:,4]

	return

#=========================================================================================
# core functions
#=========================================================================================

def calculate_avg():													#DONE

	global avg_clustering_upper_POPC
	global avg_clustering_upper_POPE
	global avg_clustering_upper_POPS
	global avg_clustering_lower_POPC
	global avg_clustering_lower_POPE
	global avg_clustering_lower_POPS
	global std_clustering_upper_POPC
	global std_clustering_upper_POPE
	global std_clustering_upper_POPS
	global std_clustering_lower_POPC
	global std_clustering_lower_POPE
	global std_clustering_lower_POPS
				
	avg_clustering_upper_POPC = np.zeros((nb_rows,1))
	avg_clustering_upper_POPE = np.zeros((nb_rows,1))
	avg_clustering_upper_POPS = np.zeros((nb_rows,1))
	avg_clustering_lower_POPC = np.zeros((nb_rows,1))
	avg_clustering_lower_POPE = np.zeros((nb_rows,1))
	avg_clustering_lower_POPS = np.zeros((nb_rows,1))
	std_clustering_upper_POPC = np.zeros((nb_rows,1))
	std_clustering_upper_POPE = np.zeros((nb_rows,1))
	std_clustering_upper_POPS = np.zeros((nb_rows,1))
	std_clustering_lower_POPC = np.zeros((nb_rows,1))
	std_clustering_lower_POPE = np.zeros((nb_rows,1))
	std_clustering_lower_POPS = np.zeros((nb_rows,1))

	#remove nan values of the weights for upper species
	#--------------------------------------------------
	#POPC
	weights_upper_POPC_nan = np.zeros((nb_rows, 1))	
	weights_upper_POPC_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_upper_POPC = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_clustering_upper_POPC[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_upper_POPC[r,0] -= 1
	weights_upper_POPC_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_upper_POPC_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_upper_POPC_nan[weights_upper_POPC_nan == 0] = 1
	
	#POPE
	weights_upper_POPE_nan = np.zeros((nb_rows, 1))	
	weights_upper_POPE_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_upper_POPE = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_clustering_upper_POPE[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_upper_POPE[r,0] -= 1	
	weights_upper_POPE_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_upper_POPE_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_upper_POPE_nan[weights_upper_POPE_nan == 0] = 1

	#POPS
	weights_upper_POPS_nan = np.zeros((nb_rows, 1))	
	weights_upper_POPS_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_upper_POPS = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_clustering_upper_POPS[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_upper_POPS[r,0] -= 1
	weights_upper_POPS_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_upper_POPS_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_upper_POPS_nan[weights_upper_POPS_nan == 0] = 1

	#remove nan values of the weights for lower species
	#--------------------------------------------------
	#POPC
	weights_lower_POPC_nan = np.zeros((nb_rows, 1))	
	weights_lower_POPC_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_lower_POPC = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_clustering_lower_POPC[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_lower_POPC[r,0] -= 1
	weights_lower_POPC_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_lower_POPC_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_lower_POPC_nan[weights_lower_POPC_nan == 0] = 1
	
	#POPE
	weights_lower_POPE_nan = np.zeros((nb_rows, 1))	
	weights_lower_POPE_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_lower_POPE = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_clustering_lower_POPE[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_lower_POPE[r,0] -= 1
	weights_lower_POPE_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_lower_POPE_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_lower_POPE_nan[weights_lower_POPE_nan == 0] = 1

	#POPS
	weights_lower_POPS_nan = np.zeros((nb_rows, 1))	
	weights_lower_POPS_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_lower_POPS = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_clustering_lower_POPS[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_lower_POPS[r,0] -= 1
	weights_lower_POPS_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_lower_POPS_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_lower_POPS_nan[weights_lower_POPS_nan == 0] = 1

	#calculate weighted average taking into account "nan"
	#----------------------------------------------------
	if args.membrane in ["AM_zCter","AM_zNter"]:
		avg_clustering_upper_POPC[:,0] =  scipy.stats.nanmean(data_clustering_upper_POPC * weights * nb_files_upper_POPC / weights_upper_POPC_nan, axis = 1)
		avg_clustering_upper_POPE[:,0] =  scipy.stats.nanmean(data_clustering_upper_POPE * weights * nb_files_upper_POPE / weights_upper_POPE_nan, axis = 1)
		avg_clustering_lower_POPC[:,0] =  scipy.stats.nanmean(data_clustering_lower_POPC * weights * nb_files_lower_POPC / weights_lower_POPC_nan, axis = 1)
		avg_clustering_lower_POPE[:,0] =  scipy.stats.nanmean(data_clustering_lower_POPE * weights * nb_files_lower_POPE / weights_lower_POPE_nan, axis = 1)
		avg_clustering_lower_POPS[:,0] =  scipy.stats.nanmean(data_clustering_lower_POPS * weights * nb_files_lower_POPS / weights_lower_POPS_nan, axis = 1)
	elif args.membrane == "SMa":
		avg_clustering_upper_POPC[:,0] =  scipy.stats.nanmean(data_clustering_upper_POPC * weights * nb_files_upper_POPC / weights_upper_POPC_nan, axis = 1)
		avg_clustering_upper_POPE[:,0] =  scipy.stats.nanmean(data_clustering_upper_POPE * weights * nb_files_upper_POPE / weights_upper_POPE_nan, axis = 1)
		avg_clustering_upper_POPS[:,0] =  scipy.stats.nanmean(data_clustering_upper_POPS * weights * nb_files_upper_POPS / weights_upper_POPS_nan, axis = 1)
		avg_clustering_lower_POPC[:,0] =  scipy.stats.nanmean(data_clustering_lower_POPC * weights * nb_files_lower_POPC / weights_lower_POPC_nan, axis = 1)
		avg_clustering_lower_POPE[:,0] =  scipy.stats.nanmean(data_clustering_lower_POPE * weights * nb_files_lower_POPE / weights_lower_POPE_nan, axis = 1)
		avg_clustering_lower_POPS[:,0] =  scipy.stats.nanmean(data_clustering_lower_POPS * weights * nb_files_lower_POPS / weights_lower_POPS_nan, axis = 1)
	elif args.membrane == "SMz":
		avg_clustering_upper_POPC[:,0] =  scipy.stats.nanmean(data_clustering_upper_POPC * weights * nb_files_upper_POPC / weights_upper_POPC_nan, axis = 1)
		avg_clustering_upper_POPE[:,0] =  scipy.stats.nanmean(data_clustering_upper_POPE * weights * nb_files_upper_POPE / weights_upper_POPE_nan, axis = 1)
		avg_clustering_lower_POPC[:,0] =  scipy.stats.nanmean(data_clustering_lower_POPC * weights * nb_files_lower_POPC / weights_lower_POPC_nan, axis = 1)
		avg_clustering_lower_POPE[:,0] =  scipy.stats.nanmean(data_clustering_lower_POPE * weights * nb_files_lower_POPE / weights_lower_POPE_nan, axis = 1)
	
	#calculate unbiased weighted std dev taking into account "nan"
	#-------------------------------------------------------------
	if args.membrane in ["AM_zCter","AM_zNter"]:
		tmp_upper_POPC = np.zeros((nb_rows, 1))
		tmp_upper_POPC[:,0] = np.nansum(weights * (data_clustering_upper_POPC - avg_clustering_upper_POPC[:,0:1])**2, axis = 1)			
		tmp_div = np.copy((weights_upper_POPC_nan)**2 - weights_upper_POPC_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_upper_POPC = np.sqrt(weights_upper_POPC_nan / tmp_div * tmp_upper_POPC)

		tmp_upper_POPE = np.zeros((nb_rows, 1))
		tmp_upper_POPE[:,0] = np.nansum(weights * (data_clustering_upper_POPE - avg_clustering_upper_POPE[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_upper_POPE_nan)**2 - weights_upper_POPE_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_upper_POPE = np.sqrt(weights_upper_POPE_nan / tmp_div * tmp_upper_POPE)

		tmp_lower_POPC = np.zeros((nb_rows, 1))
		tmp_lower_POPC[:,0] = np.nansum(weights * (data_clustering_lower_POPC - avg_clustering_lower_POPC[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_lower_POPC_nan)**2 - weights_lower_POPC_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_lower_POPC = np.sqrt(weights_lower_POPC_nan / tmp_div * tmp_lower_POPC)

		tmp_lower_POPE = np.zeros((nb_rows, 1))
		tmp_lower_POPE[:,0] = np.nansum(weights * (data_clustering_lower_POPE - avg_clustering_lower_POPE[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_lower_POPE_nan)**2 - weights_lower_POPE_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_lower_POPE = np.sqrt(weights_lower_POPE_nan / tmp_div * tmp_lower_POPE)

		tmp_lower_POPS = np.zeros((nb_rows, 1))
		tmp_lower_POPS[:,0] = np.nansum(weights * (data_clustering_lower_POPS - avg_clustering_lower_POPS[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_lower_POPS_nan)**2 - weights_lower_POPS_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_lower_POPS = np.sqrt(weights_lower_POPS_nan / tmp_div * tmp_lower_POPS)
	
	elif args.membrane == "SMa":
		tmp_upper_POPC = np.zeros((nb_rows, 1))
		tmp_upper_POPC[:,0] = np.nansum(weights * (data_clustering_upper_POPC - avg_clustering_upper_POPC[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_upper_POPC_nan)**2 - weights_upper_POPC_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_upper_POPC = np.sqrt(weights_upper_POPC_nan / tmp_div * tmp_upper_POPC)
		
		tmp_upper_POPE = np.zeros((nb_rows, 1))
		tmp_upper_POPE[:,0] = np.nansum(weights * (data_clustering_upper_POPE - avg_clustering_upper_POPE[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_upper_POPE_nan)**2 - weights_upper_POPE_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_upper_POPE = np.sqrt(weights_upper_POPE_nan / tmp_div * tmp_upper_POPE)
		
		tmp_upper_POPS = np.zeros((nb_rows, 1))
		tmp_upper_POPS[:,0] = np.nansum(weights * (data_clustering_upper_POPS - avg_clustering_upper_POPS[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_upper_POPS_nan)**2 - weights_upper_POPS_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_upper_POPS = np.sqrt(weights_upper_POPS_nan / tmp_div * tmp_upper_POPS)
		
		tmp_lower_POPC = np.zeros((nb_rows, 1))
		tmp_lower_POPC[:,0] = np.nansum(weights * (data_clustering_lower_POPC - avg_clustering_lower_POPC[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_lower_POPC_nan)**2 - weights_lower_POPC_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_lower_POPC = np.sqrt(weights_lower_POPC_nan / tmp_div * tmp_lower_POPC)
		
		tmp_lower_POPE = np.zeros((nb_rows, 1))
		tmp_lower_POPE[:,0] = np.nansum(weights * (data_clustering_lower_POPE - avg_clustering_lower_POPE[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_lower_POPE_nan)**2 - weights_lower_POPE_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_lower_POPE = np.sqrt(weights_lower_POPE_nan / tmp_div * tmp_lower_POPE)
		
		tmp_lower_POPS = np.zeros((nb_rows, 1))
		tmp_lower_POPS[:,0] = np.nansum(weights * (data_clustering_lower_POPS - avg_clustering_lower_POPS[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_lower_POPS_nan)**2 - weights_lower_POPS_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_lower_POPS = np.sqrt(weights_lower_POPS_nan / tmp_div * tmp_lower_POPS)
	
	elif args.membrane == "SMz":
		tmp_upper_POPC = np.zeros((nb_rows, 1))
		tmp_upper_POPC[:,0] = np.nansum(weights * (data_clustering_upper_POPC - avg_clustering_upper_POPC[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_upper_POPC_nan)**2 - weights_upper_POPC_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_upper_POPC = np.sqrt(weights_upper_POPC_nan / tmp_div * tmp_upper_POPC)
		
		tmp_upper_POPE = np.zeros((nb_rows, 1))
		tmp_upper_POPE[:,0] = np.nansum(weights * (data_clustering_upper_POPE - avg_clustering_upper_POPE[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_upper_POPE_nan)**2 - weights_upper_POPE_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_upper_POPE = np.sqrt(weights_upper_POPE_nan / tmp_div * tmp_upper_POPE)
		
		tmp_lower_POPC = np.zeros((nb_rows, 1))
		tmp_lower_POPC[:,0] = np.nansum(weights * (data_clustering_lower_POPC - avg_clustering_lower_POPC[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_lower_POPC_nan)**2 - weights_lower_POPC_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_lower_POPC = np.sqrt(weights_lower_POPC_nan / tmp_div * tmp_lower_POPC)

		tmp_lower_POPE = np.zeros((nb_rows, 1))
		tmp_lower_POPE[:,0] = np.nansum(weights * (data_clustering_lower_POPE - avg_clustering_lower_POPE[:,0:1])**2, axis = 1)	
		tmp_div = np.copy((weights_lower_POPE_nan)**2 - weights_lower_POPE_nan_sq)
		tmp_div[tmp_div == 0] = 1
		std_clustering_lower_POPE = np.sqrt(weights_lower_POPE_nan / tmp_div * tmp_lower_POPE)

	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():														#DONE

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [average xvg - written by xvg_average_op_simple v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - files: " + str(tmp_files[1:]) + "\n")
	if np.sum(weights) > len(args.xvgfilenames):
		output_xvg.write("# -> weight = " + str(np.sum(weights)) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Average xvg\"\n")
	output_xvg.write("@ xaxis label \"distance from cluster z axis (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"order parameter\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 12\n")
	output_xvg.write("@ s0 legend \"upper POPC (avg)\"\n")
	output_xvg.write("@ s1 legend \"upper POPC (std)\"\n")
	output_xvg.write("@ s2 legend \"upper POPE (avg)\"\n")
	output_xvg.write("@ s3 legend \"upper POPE (std)\"\n")
	output_xvg.write("@ s4 legend \"upper POPS (avg)\"\n")
	output_xvg.write("@ s5 legend \"upper POPS (std)\"\n")
	output_xvg.write("@ s6 legend \"lower POPC (avg)\"\n")
	output_xvg.write("@ s7 legend \"lower POPC (std)\"\n")
	output_xvg.write("@ s8 legend \"lower POPE (avg)\"\n")
	output_xvg.write("@ s9 legend \"lower POPE (std)\"\n")
	output_xvg.write("@ s10 legend \"lower POPS (avg)\"\n")
	output_xvg.write("@ s11 legend \"lower POPS (std)\"\n")
	
	#data
	for r in range(0, nb_rows):
		results = str(distance[r,0])
		results += "	" + "{:.6e}".format(avg_clustering_upper_POPC[r,0]) + "	" + "{:.6e}".format(std_clustering_upper_POPC[r,0]) + "	" + "{:.6e}".format(avg_clustering_upper_POPE[r,0]) + "	" + "{:.6e}".format(std_clustering_upper_POPE[r,0]) + "	" + "{:.6e}".format(avg_clustering_upper_POPS[r,0]) + "	" + "{:.6e}".format(std_clustering_upper_POPS[r,0]) + "	" + "{:.6e}".format(avg_clustering_lower_POPC[r,0]) + "	" + "{:.6e}".format(std_clustering_lower_POPC[r,0]) + "	" + "{:.6e}".format(avg_clustering_lower_POPE[r,0]) + "	" + "{:.6e}".format(std_clustering_lower_POPE[r,0]) + "	" + "{:.6e}".format(avg_clustering_lower_POPS[r,0]) + "	" + "{:.6e}".format(std_clustering_lower_POPS[r,0])
		output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

##########################################################################################
# MAIN
##########################################################################################

print "\nReading files..."
load_xvg()

print "\n\nWriting average file..."
calculate_avg()
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + args.output_file + ".xvg'."
print ""
sys.exit(0)
