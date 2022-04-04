# Loop for delta a
#
# Author: Joseph A'Hearn
# Created 01/11/2021
#
# This program runs simulations 
#   in order to make a plot of delta_a/a vs. mass


import numpy as np 
from subprocess import call
import sys
import os
import add_bodies as ab 

def load_file_list():
	file_list = ['simulation_summary.png', 'element.out', 'body0001.aei', 'body0002.aei', 'info.out', 'restart.dmp', 'param.dmp', 'small.dmp', 'big.dmp', 'restart.tmp', 'param.tmp', 'small.tmp', 'big.tmp', 'xv.out', 'ce.out']
	return file_list

def plot_results():
	# needs to be implemented

def main(min_mass_exp=-19, max_mass_exp=-10, min_delta_a=0.01, max_delta_a=0.1, delta_a_increment=0.01): # delta a in units of km 
#def main(min_mass_exp=-19, max_mass_exp=-10, min_delta_a=0.01, max_delta_a=1.0, delta_a_increment=0.05): # delta a in units of km 
	path = '~/Downloads/03.Coorbital_Projects/D68_model/'
	file_list = load_file_list()
	call('rm *mp *out *aei', shell=True)

	for j in range(min_mass_exp, max_mass_exp, 2):
		for i in range(min_delta_a, max_delta_a, delta_a_increment): # i corresponds to the a_perturbation
			ab.Obj5_leading_D68_configuration(a_perturbation=min_delta_a)
			call('./mercury6', shell=True)
			print('About to run ./element6')
			call('./element6', shell=True)

			print('About to run simsum.py') # still need to modify simsum.py to make this work 
			call('pythonw simsum.py', shell=True)

			# save plot and other files in a new folder
			folder_name = 'simulation_' + str(j) + '_' + str(i) 
			call('mkdir ' + folder_name, shell=True)

			#call('rm a_vs_t.png', shell=True)

			for file in file_list:
				call('mv ' + path + file + ' ' + path + folder_name + '/' + file, shell=True)
		print('Mass 10^' + str(j) + ' completed.')
	print('Entire loop completed.')

	plot_results()