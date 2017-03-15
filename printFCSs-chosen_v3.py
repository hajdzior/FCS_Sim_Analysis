import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pylab import *
import os
from os.path import isfile, join
import re # Regular expres.

import pdb



#MW scan folders and subfolders looking for paths containing: something-number-dot-dat and save the path
def scan_dir(dir):
	list_dat = []
	for name in os.listdir(dir):
		path = os.path.join(dir, name)
		if os.path.isfile(path):
			#if re.search('.+\d+\.dat', path): # $:end-line mark
				#path_name = re.search('(.+\d+)\.dat$', path).group(0)
			if re.search('.*.ACF.dat', path): # $:end-line mark
				path_name = re.search('.*.ACF.dat$', path).group(0)
				list_dat.append(path_name)
		else:
			list_dat = list_dat + scan_dir(path)
	return list_dat

#def find_logFile(dir):
	#list_dat = []
	#for name in os.listdir(dir):
		#path = os.path.join(dir, name)
		#if os.path.isfile(path):
			#if re.search('\d+\.dat$', path): # $:end-line mark
				#path_name = re.search('(.+\d+)\.dat$', path).group(0)
				#list_dat.append(path_name)
		#else:
			#list_dat = list_dat + find_logFile(path)
	#return list_dat

main_dir = "./SIM_data"
list_dir = os.listdir(main_dir)
print list_dir

k_ex_mxs = []
for dir_now in list_dir:
	print dir_now
	
	# Read FSC data
	print 'Im reading FCS data'
	scan_now = main_dir + '/' + dir_now
	paths_to_data = scan_dir(scan_now)
	print paths_to_data 
	
	# set colors
	cmap = plt.get_cmap('jet_r')
	N = len(paths_to_data)
	c_num = -1.
	fig = figure(1, figsize=(10,7)) #
	
	#for data_now in range(len(paths_to_data) - 1):	
	for data_now in range(0, len(paths_to_data)):	
		inputfile = paths_to_data[data_now]
		data = np.loadtxt(inputfile, skiprows=0)
		
		#read log file
		log_ind = inputfile.find('/data/') 
		log_file = inputfile[0:log_ind] + '/log.txt'
		f = open (log_file, "rt")
		l = f.readline()
		while l:
			va = l.split()
			if (va[0] == 'k_ex_max='):
				k_ex_mxs = float(va[-1])
			l = f.readline()
		f.close()
		print k_ex_mxs
		#print inputfile
		# colors
		c_num = c_num + 1.
		color = cmap(c_num/N)
		fig=figure(1)
		lab = k_ex_mxs
		FCS_data = data.transpose()
		
		#plot FCS func
		plt.plot(FCS_data[0], FCS_data[1], label = lab, c = color)

		
	# plot settings
	plt.legend(loc='best')
	#plt.xlim(0.001,100.)
	plt.xscale('log')
	plt.ylim(-0.,2.2)
	plt.xlim(0.0000001,0.1)
	plt.ylabel('FCS')
	plt.xlabel('$\\tau$ [s]')
	plt.legend(ncol=1, loc='best', 
	columnspacing=1.0, labelspacing=0.0,
	handletextpad=0.0, handlelength=1.5,
	fancybox=True, shadow=True)# bbox_to_anchor=[0.5, 1.1]
	f_name = 'FCSs/' + dir_now + '.png'
	plt.savefig(f_name)
	plt.clf()

	

 
