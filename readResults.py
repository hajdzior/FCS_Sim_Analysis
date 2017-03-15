#!/usr/bin/python
# -*- coding: utf-8 -*-
"""  
"""
#from __future__ import print_function

import numpy as np
import os
from os.path import abspath, dirname, join
import sys
import time
import re

from matplotlib import pylab as plt

#exec(open("./config.cfg").read(), globals()) 

#sys.path.insert(0, dirname(dirname(abspath(__file__))))

from multipletau import autocorrelate #, correlate, correlate_numpy

#MW scan folders and subfolders looking for paths containing: something-number-dot-dat and save the path
def scan_dir(dir):
	list_dat = []
	for name in os.listdir(dir):
		path = os.path.join(dir, name)
		if os.path.isfile(path):
			if re.search('\d+\.dat$', path): # $:end-line mark
				path_name = re.search('(.+\d+)\.dat$', path).group(0)
				list_dat.append(path_name)
		else:
			list_dat = list_dat + scan_dir(path)
	return list_dat


main_dir = "./SIM_data"
list_dir = os.listdir(main_dir)
print list_dir

#colors
cmap = plt.get_cmap('jet_r')
N = len(list_dir)
c_num = -1.

plt.figure(1, figsize=(10,7))
for dir_now in list_dir:
	print dir_now
	paths_to_data = scan_dir(main_dir + '/' + dir_now)

	k_ex_mxs = []
	kcpss = []
	print paths_to_data
	
	#set color for this data set
	c_num = c_num + 1.
	color = cmap(c_num/N)


	for file_name in paths_to_data:
		log_ind = file_name.find('/data/') 
		log_file = file_name[0:log_ind] + '/log.txt'
		print log_file
		
		f = open (log_file, "rt")
		l = f.readline()
		while l:
			va = l.split()
			if (va[0] == 'k_ex_max='):
				k_ex_mxs.append(float(va[-1]))
			if (va[0] == 'avg'):
				kcpss.append(float(va[3]))
			l = f.readline()
		f.close()
	print k_ex_mxs
	print kcpss
	
	lab_now = dir_now
	
	plt.plot(k_ex_mxs, kcpss, linestyle="", marker="o", color=color, label = lab_now)
		
	plt.xlabel('k_ex_max', fontsize = 20)
	plt.ylabel('kcps', fontsize = 20)
	plt.rc('xtick',labelsize=20) 
	plt.rc('ytick',labelsize=20) 
	plt.legend(loc = 'best', fontsize=14)
	f_name = './ResPlots/intVSexc.png'
	plt.savefig(f_name)

