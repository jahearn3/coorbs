# Overarching Plots
#
# Author: Joseph A'Hearn
# Created 03/11/2020
#
# This program plots data from many simulations
#   

import numpy as np 
import plot_assistant as pa 
import matplotlib.pyplot as plt 
import bodies as bd 

# =================================================================================================================

# IMPORTING FUNCTIONS

def import_data_2():
	data = np.loadtxt('stats_2.txt', skiprows=1)
	M            = data[:,0] 						 
	delta_a      = data[:,1] 
	delta_lambda = data[:,2] 
	min_sep      = data[:,3] 
	initial_sep  = data[:,4] 
	return M, delta_a, delta_lambda, min_sep, initial_sep

def import_data_3():
	data = np.loadtxt('stats_3.txt', skiprows=1)
	M            = data[:,0] 						 
	delta_a      = data[:,1] 
	min_sep1     = data[:,2] 
	initial_sep1 = data[:,3] 
	min_sep2     = data[:,4] 
	initial_sep2 = data[:,5] 
	dlmda1       = data[:,6] 
	dlmda2       = data[:,7] 
	dlmda3       = data[:,8] 
	return M, delta_a, min_sep1, initial_sep1, min_sep2, initial_sep2, dlmda1, dlmda2, dlmda3 

def import_data_4():
	data = np.loadtxt('stats_4.txt', skiprows=1)
	M            = data[:,0] 						 
	delta_a      = data[:,1] 
	delta_lambda = data[:,2] 
	min_sep1     = data[:,3] 
	initial_sep1 = data[:,4] 
	min_sep2     = data[:,5] 
	initial_sep2 = data[:,6] 
	min_sep3     = data[:,7] 
	initial_sep3 = data[:,8] 
	return M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3

def import_data_5():
	data = np.loadtxt('stats_5.txt', skiprows=1)
	M            = data[:,0] 						 
	delta_a      = data[:,1] 
	delta_lambda = data[:,2] 
	min_sep1     = data[:,3] 
	initial_sep1 = data[:,4] 
	min_sep2     = data[:,5] 
	initial_sep2 = data[:,6] 
	min_sep3     = data[:,7] 
	initial_sep3 = data[:,8] 
	min_sep4     = data[:,9] 
	initial_sep4 = data[:,10] 
	return M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3, min_sep4, initial_sep4

# =================================================================================================================

# PLOTTING FUNCTIONS

# N = 2

def plot_delta_a_by_initial_separation_2(iniseps, radialoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'radial oscillation amplitude [km]', 0, 185)
	for i in range(3, len(iniseps)):
		ax.plot(iniseps[i], radialoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.set_xticks([0, 60, 120, 180])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=1)
	ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='radial_oscillation_by_initial_sep_2', black=black)

def plot_delta_a_by_mass_2(iniseps, radialoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'radial oscillation amplitude [km]')
	for i in range(len(masses)):
		ax.plot(masses[i], radialoscs[i], label=str(iniseps[i][0]))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=0, fontsize=5)
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='radial_oscillation_by_mass_2', black=black)

def plot_delta_lambda_by_initial_separation_2(iniseps, longitoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'longitudinal oscillation [deg]', 0, 185)
	for i in range(3, len(iniseps)):
		ax.plot(iniseps[i], longitoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.set_xticks([0, 60, 120, 180])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=1)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitudinal_oscillation_by_initial_sep_2', black=black)

def plot_delta_lambda_by_mass_2(iniseps, longitoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'longitudinal oscillation [deg]')
	for i in range(len(iniseps)):
		ax.plot(masses[i], longitoscs[i], label=str(iniseps[i][0]))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.set_xscale('log')
	ax.legend(loc=0, fontsize=5)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitudinal_oscillation_by_mass_2', black=black)

def plot_minimum_separation_by_initial_separation_2(iniseps, minseps, masses, black=False): 
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'minimum separation [deg]', 0, 185, 0, 185)
	#ax.scatter(initial_sep, min_sep)
	#for i in range(len(iniseps)):
		#ax.scatter(iniseps[i], minseps[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.scatter(iniseps[0], minseps[0], label=r'$\geq$' + "{:.1e}".format(masses[0][0]) + ' kg')
	#ax.scatter(iniseps[1], minseps[1], label="{:.1e}".format(masses[1][0]) + ' kg')
	#ax.scatter(iniseps[2], minseps[2], label="{:.1e}".format(masses[2][0]) + ' kg')
	#ax.scatter(iniseps[3], minseps[3], label="{:.1e}".format(masses[3][0]) + ' kg')
	ax.scatter(iniseps[4], minseps[4], label=r'$\leq$' + "{:.1e}".format(masses[4][0]) + ' kg') 
	ax.scatter(iniseps[5], minseps[5], label="{:.1e}".format(masses[5][0]) + ' kg') # middle line
	#ax.scatter(iniseps[6], minseps[6], label="{:.1e}".format(masses[6][0]) + ' kg')
	#ax.scatter(iniseps[7], minseps[7], label="{:.1e}".format(masses[7][0]) + ' kg')
	#ax.scatter(iniseps[8], minseps[8], label="{:.1e}".format(masses[8][0]) + ' kg')
	ax.set_xticks([0, 60, 120, 180])
	ax.set_yticks([0, 60, 120, 180])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=2)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='min_sep_by_initial_sep_2', black=black)

def plot_minimum_separation_by_mass_2(iniseps, minseps, masses, black=False): 
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'minimum separation [deg]')
	#ax.scatter(initial_sep, min_sep)
	for i in range(len(iniseps)):
		ax.plot(masses[i], minseps[i], label=str(iniseps[i][0]))
	ax.set_xscale('log')
	ax.set_yticks([0, 60, 120, 180])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=0, fontsize=5)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='min_sep_by_mass_2', black=black)

# ============================================================================================================

# N = 3

def plot_delta_a_by_initial_separation_3(inisep1s, inisep2s, radialoscs, masses, M=bd.M(), black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'radial oscillation amplitude [' + r'$R_H$' + ']')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(inisep1s[i], radialoscs[i], label='mass ratio ' + "{:.1e}".format(masses[i][0] / M))
		# inisep2 must be equal to inisep1 in a symmetric configuration, so no need to try to plot it again
		#ax.plot(inisep2s[i], radialoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.set_xticks([0, 60, 120, 180])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=4)
	#ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='radial_oscillation_by_initial_sep_3', black=black)

def plot_delta_a_by_mass_3(inisep1s, inisep2s, radialoscs, masses, M=bd.M(), black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'radial oscillation amplitude [' + r'$R_H$' + ']')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(masses))))

	for i in range(len(masses)):
		ax.scatter(masses[i][0] / M, radialoscs[i], label=str(float(inisep1s[i][0])))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=0, fontsize=5, ncol=4)
	ax.set_xscale('log')
	#ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='radial_oscillation_by_mass_3', black=black)

def plot_delta_lambda_by_initial_separation_3(inisep1s, inisep2s, longitosc1s, longitosc2s, longitosc3s, masses, M=bd.M(), black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'longitudinal oscillation [deg]')
	#ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.scatter(inisep1s[i], longitosc1s[i], label='mass ratio ' + "{:.1e}".format(masses[i][0] / M))
		ax.scatter(inisep1s[i], longitosc2s[i], label='mass ratio ' + "{:.1e}".format(masses[i][0] / M))
		ax.scatter(inisep1s[i], longitosc3s[i], label='mass ratio ' + "{:.1e}".format(masses[i][0] / M))
		# inisep2 must be equal to inisep1 in a symmetric configuration, so no need to try to plot it again
		#ax.plot(inisep2s[i], longitoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.set_xticks([0, 60, 120, 180])
	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=2)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitudinal_oscillation_by_initial_sep_3', black=black)

def plot_delta_lambda_by_mass_3(inisep1s, inisep2s, longitosc1s, longitosc2s, longitosc3s, masses, M=bd.M(), black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'longitudinal oscillation [deg]')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.scatter(masses[i] / M, longitosc1s[i], label=str(inisep1s[i][0]))
		ax.scatter(masses[i] / M, longitosc2s[i], label=str(inisep1s[i][0]))
		ax.scatter(masses[i] / M, longitosc3s[i], label=str(inisep1s[i][0]))
		#ax.plot(masses[i], longitoscs[i], label=str(inisep2s[i][0]))
	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.set_xscale('log')
	ax.legend(loc=0, fontsize=5, ncol=5)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitudinal_oscillation_by_mass_3', black=black)

def plot_minimum_separation_by_initial_separation_3(inisep1s, minsep1s, inisep2s, minsep2s, masses, M=bd.M(), black=False): 
	minminseps = np.minimum(minsep1s, minsep2s)
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'minimum separation [deg]')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(inisep1s[i], minminseps[i], label='mass ratio ' + "{:.1e}".format(masses[i][0] / M))
		#ax.plot(inisep2s[i], minsep2s[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	#ax.scatter(iniseps[0], minseps[0], label=r'$\geq$' + "{:.1e}".format(masses[0][0]) + ' kg')
	#ax.scatter(iniseps[1], minseps[1], label="{:.1e}".format(masses[1][0]) + ' kg')
	#ax.scatter(iniseps[2], minseps[2], label="{:.1e}".format(masses[2][0]) + ' kg')
	#ax.scatter(iniseps[3], minseps[3], label="{:.1e}".format(masses[3][0]) + ' kg')
	#ax.scatter(iniseps[4], minseps[4], label=r'$\leq$' + "{:.1e}".format(masses[4][0]) + ' kg') 
	#ax.scatter(iniseps[5], minseps[5], label="{:.1e}".format(masses[5][0]) + ' kg') # middle line
	#ax.scatter(iniseps[6], minseps[6], label="{:.1e}".format(masses[6][0]) + ' kg')
	#ax.scatter(iniseps[7], minseps[7], label="{:.1e}".format(masses[7][0]) + ' kg')
	#ax.scatter(iniseps[8], minseps[8], label="{:.1e}".format(masses[8][0]) + ' kg')
	ax.set_xticks([0, 60, 120, 180])
	ax.set_yticks([0, 30, 60, 90, 120])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=2)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='min_sep_by_initial_sep_3', black=black)

def plot_minimum_separation_by_mass_3(inisep1s, minsep1s, inisep2s, minsep2s, masses, M=bd.M(), black=False): 
	minminseps = np.minimum(minsep1s, minsep2s)
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'minimum separation [deg]')
	#ax.scatter(initial_sep, min_sep)
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.scatter(masses[i] / M, minminseps[i], label=str(inisep1s[i][0]))
		#ax.plot(masses[i], minsep2s[i], label=str(inisep2s[i][0]))
	#ax.scatter(iniseps[0], minseps[0], label=r'$\geq$' + "{:.1e}".format(masses[0][0]) + ' kg')
	#ax.scatter(iniseps[1], minseps[1], label="{:.1e}".format(masses[1][0]) + ' kg')
	#ax.scatter(iniseps[2], minseps[2], label="{:.1e}".format(masses[2][0]) + ' kg')
	#ax.scatter(iniseps[3], minseps[3], label="{:.1e}".format(masses[3][0]) + ' kg')
	#ax.scatter(iniseps[4], minseps[4], label=r'$\leq$' + "{:.1e}".format(masses[4][0]) + ' kg') 
	#ax.scatter(iniseps[5], minseps[5], label="{:.1e}".format(masses[5][0]) + ' kg') # middle line
	#ax.scatter(iniseps[6], minseps[6], label="{:.1e}".format(masses[6][0]) + ' kg')
	#ax.scatter(iniseps[7], minseps[7], label="{:.1e}".format(masses[7][0]) + ' kg')
	#ax.scatter(iniseps[8], minseps[8], label="{:.1e}".format(masses[8][0]) + ' kg')
	ax.set_xscale('log')
	ax.set_yticks([0, 60, 120, 180])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=0, fontsize=5, ncol=6)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='min_sep_by_mass_3', black=black)

# =================================================================================================================

# N = 4

def plot_delta_a_by_initial_separation_4(inisep1s, inisep2s, inisep3s, radialoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'radial oscillation amplitude [km]')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(inisep1s[i], radialoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
		#ax.plot(inisep2s[i], radialoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.set_xticks([20, 40, 60, 80, 100])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=4)
	ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='radial_oscillation_by_initial_sep_4', black=black)

def plot_delta_a_by_mass_4(inisep1s, inisep2s, inisep3s, radialoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'radial oscillation amplitude [km]')
	#print(len(masses))
	#print(len(inisep1s))
	#print(len(inisep1s[0]))
	#print(len(masses[0]))
	#print(inisep1s[0][0])
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(masses))))
	for i in range(len(masses)):
		#print(len(masses[i]))
		#print(len(radialoscs[i]))
		#print(len(inisep1s[i]))
		#print(str(float(inisep1s[i][0])))
		ax.plot(masses[i], radialoscs[i], label=str(float(inisep1s[i][0])))
		#ax.plot(masses[i], radialoscs[i], label=str(float(inisep2s[i][0])))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=0, fontsize=5, ncol=4)
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='radial_oscillation_by_mass_4', black=black)

def plot_delta_lambda_by_initial_separation_4(inisep1s, inisep2s, inisep3s, longitoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'longitudinal oscillation [deg]')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(inisep1s[i], longitoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
		#ax.plot(inisep2s[i], longitoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.set_xticks([20, 40, 60, 80, 100])
	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=3)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitudinal_oscillation_by_initial_sep_4', black=black)

def plot_delta_lambda_by_mass_4(inisep1s, inisep2s, inisep3s, longitoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'longitudinal oscillation [deg]')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(masses[i], longitoscs[i], label=str(inisep1s[i][0]))
		#ax.plot(masses[i], longitoscs[i], label=str(inisep2s[i][0]))
	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.set_xscale('log')
	ax.legend(loc=0, fontsize=5, ncol=5)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitudinal_oscillation_by_mass_4', black=black)

def plot_minimum_separation_by_initial_separation_4(inisep1s, minsep1s, inisep2s, minsep2s, inisep3s, minsep3s, masses, black=False): 
	min12seps  = np.minimum(minsep1s,  minsep2s)
	minminseps = np.minimum(min12seps, minsep3s)
	# one issue here is that the initial separation is not necessarily equal for all bodies
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'minimum separation [deg]')
	#ax.scatter(initial_sep, min_sep)
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(inisep1s[i], minminseps[i], label="{:.1e}".format(masses[i][0]) + ' kg')
		#ax.plot(inisep2s[i], minsep2s[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	#ax.scatter(iniseps[0], minseps[0], label=r'$\geq$' + "{:.1e}".format(masses[0][0]) + ' kg')
	#ax.scatter(iniseps[1], minseps[1], label="{:.1e}".format(masses[1][0]) + ' kg')
	#ax.scatter(iniseps[2], minseps[2], label="{:.1e}".format(masses[2][0]) + ' kg')
	#ax.scatter(iniseps[3], minseps[3], label="{:.1e}".format(masses[3][0]) + ' kg')
	#ax.scatter(iniseps[4], minseps[4], label=r'$\leq$' + "{:.1e}".format(masses[4][0]) + ' kg') 
	#ax.scatter(iniseps[5], minseps[5], label="{:.1e}".format(masses[5][0]) + ' kg') # middle line
	#ax.scatter(iniseps[6], minseps[6], label="{:.1e}".format(masses[6][0]) + ' kg')
	#ax.scatter(iniseps[7], minseps[7], label="{:.1e}".format(masses[7][0]) + ' kg')
	#ax.scatter(iniseps[8], minseps[8], label="{:.1e}".format(masses[8][0]) + ' kg')
	ax.set_xticks([20, 40, 60, 80, 100])
	ax.set_yticks([0, 10, 20, 30, 40])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=3)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='min_sep_by_initial_sep_4', black=black)

def plot_minimum_separation_by_mass_4(inisep1s, minsep1s, inisep2s, minsep2s, inisep3s, minsep3s, masses, black=False): 
	min12seps  = np.minimum(minsep1s,  minsep2s)
	minminseps = np.minimum(min12seps, minsep3s)
	# one issue here is that the initial separation is not necessarily equal for all bodies
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'minimum separation [deg]')
	#ax.scatter(initial_sep, min_sep)
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(masses[i], minminseps[i], label=str(inisep1s[i][0]))
		#ax.plot(masses[i], minsep2s[i], label=str(inisep2s[i][0]))
	#ax.scatter(iniseps[0], minseps[0], label=r'$\geq$' + "{:.1e}".format(masses[0][0]) + ' kg')
	#ax.scatter(iniseps[1], minseps[1], label="{:.1e}".format(masses[1][0]) + ' kg')
	#ax.scatter(iniseps[2], minseps[2], label="{:.1e}".format(masses[2][0]) + ' kg')
	#ax.scatter(iniseps[3], minseps[3], label="{:.1e}".format(masses[3][0]) + ' kg')
	#ax.scatter(iniseps[4], minseps[4], label=r'$\leq$' + "{:.1e}".format(masses[4][0]) + ' kg') 
	#ax.scatter(iniseps[5], minseps[5], label="{:.1e}".format(masses[5][0]) + ' kg') # middle line
	#ax.scatter(iniseps[6], minseps[6], label="{:.1e}".format(masses[6][0]) + ' kg')
	#ax.scatter(iniseps[7], minseps[7], label="{:.1e}".format(masses[7][0]) + ' kg')
	#ax.scatter(iniseps[8], minseps[8], label="{:.1e}".format(masses[8][0]) + ' kg')
	ax.set_xscale('log')
	ax.set_yticks([0, 10, 20, 30, 40, 50])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=0, fontsize=5, ncol=6)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='min_sep_by_mass_4', black=black)

# =================================================================================================================

# N = 5

def plot_delta_a_by_initial_separation_5(inisep1s, inisep2s, inisep3s, inisep4s, radialoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'radial oscillation amplitude [km]')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(inisep1s[i], radialoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
		#ax.plot(inisep2s[i], radialoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.set_xticks([0, 20, 40, 60, 80, 100])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=3)
	ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='radial_oscillation_by_initial_sep_5', black=black)

def plot_delta_a_by_mass_5(inisep1s, inisep2s, inisep3s, inisep4s, radialoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'radial oscillation amplitude [km]')
	#print(len(masses))
	#print(len(inisep1s))
	#print(len(inisep1s[0]))
	#print(len(masses[0]))
	#print(inisep1s[0][0])
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(masses))))
	for i in range(len(masses)):
		#print(len(masses[i]))
		#print(len(radialoscs[i]))
		#print(len(inisep1s[i]))
		#print(str(float(inisep1s[i][0])))
		ax.plot(masses[i], radialoscs[i], label=str(float(inisep1s[i][0])))
		#ax.plot(masses[i], radialoscs[i], label=str(float(inisep2s[i][0])))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=0, fontsize=5, ncol=4)
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='radial_oscillation_by_mass_5', black=black)

def plot_delta_lambda_by_initial_separation_5(inisep1s, inisep2s, inisep3s, inisep4s, longitoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'longitudinal oscillation [deg]')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(inisep1s[i], longitoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
		#ax.plot(inisep2s[i], longitoscs[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	ax.set_xticks([0, 20, 40, 60, 80, 100])
	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=3)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitudinal_oscillation_by_initial_sep_5', black=black)

def plot_delta_lambda_by_mass_5(inisep1s, inisep2s, inisep3s, inisep4s, longitoscs, masses, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'longitudinal oscillation [deg]')
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(masses[i], longitoscs[i], label=str(inisep1s[i][0]))
		#ax.plot(masses[i], longitoscs[i], label=str(inisep2s[i][0]))
	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.set_xscale('log')
	ax.legend(loc=0, fontsize=5, ncol=5)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitudinal_oscillation_by_mass_5', black=black)

def plot_minimum_separation_by_initial_separation_5(inisep1s, minsep1s, inisep2s, minsep2s, inisep3s, minsep3s, inisep4s, minsep4s, masses, black=False): 
	min12seps = np.minimum(minsep1s, minsep2s)
	min34seps = np.minimum(minsep3s, minsep4s)
	minminseps = np.minimum(min12seps, min34seps)
	# one issue here is that the initial separation is not necessarily equal for all bodies
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'initial separation [deg]', 'minimum separation [deg]')
	#ax.scatter(initial_sep, min_sep)
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(inisep1s[i], minminseps[i], label="{:.1e}".format(masses[i][0]) + ' kg')
		#ax.plot(inisep2s[i], minsep2s[i], label="{:.1e}".format(masses[i][0]) + ' kg')
	#ax.scatter(iniseps[0], minseps[0], label=r'$\geq$' + "{:.1e}".format(masses[0][0]) + ' kg')
	#ax.scatter(iniseps[1], minseps[1], label="{:.1e}".format(masses[1][0]) + ' kg')
	#ax.scatter(iniseps[2], minseps[2], label="{:.1e}".format(masses[2][0]) + ' kg')
	#ax.scatter(iniseps[3], minseps[3], label="{:.1e}".format(masses[3][0]) + ' kg')
	#ax.scatter(iniseps[4], minseps[4], label=r'$\leq$' + "{:.1e}".format(masses[4][0]) + ' kg') 
	#ax.scatter(iniseps[5], minseps[5], label="{:.1e}".format(masses[5][0]) + ' kg') # middle line
	#ax.scatter(iniseps[6], minseps[6], label="{:.1e}".format(masses[6][0]) + ' kg')
	#ax.scatter(iniseps[7], minseps[7], label="{:.1e}".format(masses[7][0]) + ' kg')
	#ax.scatter(iniseps[8], minseps[8], label="{:.1e}".format(masses[8][0]) + ' kg')
	ax.set_xticks([0, 20, 40, 60, 80, 100])
	ax.set_yticks([0, 10, 20, 30, 40])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=2)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='min_sep_by_initial_sep_5', black=black)

def plot_minimum_separation_by_mass_5(inisep1s, minsep1s, inisep2s, minsep2s, inisep3s, minsep3s, inisep4s, minsep4s, masses, black=False): 
	min12seps = np.minimum(minsep1s, minsep2s)
	min34seps = np.minimum(minsep3s, minsep4s)
	minminseps = np.minimum(min12seps, min34seps)
	# one issue here is that the initial separation is not necessarily equal for all bodies
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'mass [kg]', 'minimum separation [deg]')
	#ax.scatter(initial_sep, min_sep)
	ax.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(inisep1s))))
	for i in range(len(inisep1s)):
		ax.plot(masses[i], minminseps[i], label=str(inisep1s[i][0]))
		#ax.plot(masses[i], minsep2s[i], label=str(inisep2s[i][0]))
	#ax.scatter(iniseps[0], minseps[0], label=r'$\geq$' + "{:.1e}".format(masses[0][0]) + ' kg')
	#ax.scatter(iniseps[1], minseps[1], label="{:.1e}".format(masses[1][0]) + ' kg')
	#ax.scatter(iniseps[2], minseps[2], label="{:.1e}".format(masses[2][0]) + ' kg')
	#ax.scatter(iniseps[3], minseps[3], label="{:.1e}".format(masses[3][0]) + ' kg')
	#ax.scatter(iniseps[4], minseps[4], label=r'$\leq$' + "{:.1e}".format(masses[4][0]) + ' kg') 
	#ax.scatter(iniseps[5], minseps[5], label="{:.1e}".format(masses[5][0]) + ' kg') # middle line
	#ax.scatter(iniseps[6], minseps[6], label="{:.1e}".format(masses[6][0]) + ' kg')
	#ax.scatter(iniseps[7], minseps[7], label="{:.1e}".format(masses[7][0]) + ' kg')
	#ax.scatter(iniseps[8], minseps[8], label="{:.1e}".format(masses[8][0]) + ' kg')
	ax.set_xscale('log')
	ax.set_yticks([0, 10, 20, 30, 40, 50])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', direction='in')
	ax.legend(loc=0, fontsize=5, ncol=6)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='min_sep_by_mass_5', black=black)

# =================================================================================================================

# CATEGORIZATION FUNCTIONS

def categorize_by_masses_2(M, delta_a, delta_lambda, min_sep, initial_sep):
	minseps = []
	iniseps = []
	masses = []
	radialoscs = []
	longitoscs = []

	minsepcur = []
	inisepcur = []
	masscur = []
	radosccur = []
	lonosccur = []


	minsepcur.append(min_sep[0])
	inisepcur.append(initial_sep[0])
	masscur.append(M[0])
	radosccur.append(delta_a[0])
	lonosccur.append(delta_lambda[0])

	for i in range(1, len(min_sep)):
		if(M[i] != M[i-1]):
			minseps.append(minsepcur)
			iniseps.append(inisepcur)
			masses.append(masscur)
			radialoscs.append(radosccur)
			longitoscs.append(lonosccur)


			minsepcur = []
			inisepcur = []
			masscur = []
			radosccur = []
			lonosccur = []


		minsepcur.append(min_sep[i])
		inisepcur.append(initial_sep[i])
		masscur.append(M[i])
		radosccur.append(delta_a[i])
		lonosccur.append(delta_lambda[i])
	else:
		minseps.append(minsepcur)
		iniseps.append(inisepcur)
		masses.append(masscur)
		radialoscs.append(radosccur)
		longitoscs.append(lonosccur)

	return iniseps, minseps, masses, radialoscs, longitoscs 

def categorize_by_inisep_2(M, delta_a, delta_lambda, min_sep, initial_sep):
	minseps = []
	iniseps = []
	masses = []
	radialoscs = []
	longitoscs = []

	for j in range(5, 185, 5):
		minsepcur = []
		inisepcur = []
		masscur = []
		radosccur = []
		lonosccur = []
		for i in range(len(min_sep)):
			if(j == initial_sep[i]):
				minsepcur.append(min_sep[i])
				inisepcur.append(initial_sep[i])
				masscur.append(M[i])
				radosccur.append(delta_a[i])
				lonosccur.append(delta_lambda[i])
		minseps.append(minsepcur)
		iniseps.append(inisepcur)
		masses.append(masscur)
		radialoscs.append(radosccur)
		longitoscs.append(lonosccur)

	return iniseps, minseps, masses, radialoscs, longitoscs

def categorize_by_masses_3(M, delta_a, dlmda1, dlmda2, dlmda3, min_sep1, initial_sep1, min_sep2, initial_sep2):
	minsep1s = []
	inisep1s = []
	minsep2s = []
	inisep2s = []
	masses = []
	radialoscs = []
	longitosc1s = []
	longitosc2s = []
	longitosc3s = []

	minsep1cur = []
	inisep1cur = []
	minsep2cur = []
	inisep2cur = []
	masscur = []
	radosccur = []
	lonosc1cur = []
	lonosc2cur = []
	lonosc3cur = []


	minsep1cur.append(min_sep1[0])
	inisep1cur.append(initial_sep1[0])
	minsep2cur.append(min_sep2[0])
	inisep2cur.append(initial_sep2[0])
	masscur.append(M[0])
	radosccur.append(delta_a[0])
	lonosc1cur.append(dlmda1[0])
	lonosc2cur.append(dlmda2[0])
	lonosc3cur.append(dlmda3[0])

	for i in range(1, len(min_sep1)):
		if(M[i] != M[i-1]):
			minsep1s.append(minsep1cur)
			inisep1s.append(inisep1cur)
			minsep2s.append(minsep2cur)
			inisep2s.append(inisep2cur)
			masses.append(masscur)
			radialoscs.append(radosccur)
			longitosc1s.append(lonosc1cur)
			longitosc2s.append(lonosc2cur)
			longitosc3s.append(lonosc3cur)


			minsep1cur = []
			inisep1cur = []
			minsep2cur = []
			inisep2cur = []
			masscur = []
			radosccur = []
			lonosc1cur = []
			lonosc2cur = []
			lonosc3cur = []


		minsep1cur.append(min_sep1[i])
		inisep1cur.append(initial_sep1[i])
		minsep2cur.append(min_sep2[i])
		inisep2cur.append(initial_sep2[i])
		masscur.append(M[i])
		radosccur.append(delta_a[i])
		lonosc1cur.append(dlmda1[i])
		lonosc2cur.append(dlmda2[i])
		lonosc3cur.append(dlmda3[i])
	else:
		minsep1s.append(minsep1cur)
		inisep1s.append(inisep1cur)
		minsep2s.append(minsep2cur)
		inisep2s.append(inisep2cur)
		masses.append(masscur)
		radialoscs.append(radosccur)
		longitosc1s.append(lonosc1cur)
		longitosc2s.append(lonosc2cur)
		longitosc3s.append(lonosc3cur)

	return inisep1s, minsep1s, inisep2s, minsep2s, masses, radialoscs, longitosc1s, longitosc2s, longitosc3s 

def categorize_by_inisep_3(M, delta_a, dlmda1, dlmda2, dlmda3, min_sep1, initial_sep1, min_sep2, initial_sep2):
	minsep1s = []
	inisep1s = []
	minsep2s = []
	inisep2s = []
	masses = []
	radialoscs = []
	longitosc1s = []
	longitosc2s = []
	longitosc3s = []

	for j in np.linspace(7.361, 179.861, 70): 
		minsep1cur = []
		inisep1cur = []
		minsep2cur = []
		inisep2cur = []
		masscur = []
		radosccur = []
		lonosc1cur = []
		lonosc2cur = []
		lonosc3cur = []
		#print(j, np.round(j, 3))
		for i in range(len(min_sep1)):
			#if(abs(j - initial_sep1[i]) < 1.0):
				#print(np.round(j, 3), np.round(initial_sep1[i], 3))
			if(np.round(j, 3) == np.round(initial_sep1[i], 3)):
				minsep1cur.append(min_sep1[i])
				inisep1cur.append(str(np.round(initial_sep1[i], 3)))
				minsep2cur.append(min_sep2[i])
				inisep2cur.append(str(np.round(initial_sep2[i], 3)))
				masscur.append(M[i])
				radosccur.append(delta_a[i])
				lonosc1cur.append(dlmda1[i])
				lonosc2cur.append(dlmda2[i])
				lonosc3cur.append(dlmda3[i])
		minsep1s.append(minsep1cur)
		inisep1s.append(inisep1cur)
		minsep2s.append(minsep2cur)
		inisep2s.append(inisep2cur)
		masses.append(masscur)
		radialoscs.append(radosccur)
		longitosc1s.append(lonosc1cur)
		longitosc2s.append(lonosc2cur)
		longitosc3s.append(lonosc3cur)

	minsep12 = []
	inisep12 = []
	minsep22 = []
	inisep22 = []
	masse2 = []
	radialosc2 = []
	longitosc12 = []
	longitosc22 = []
	longitosc32 = []

	for k in range(len(inisep1s)):
		if(len(inisep1s[k]) > 0):
			minsep12.append(minsep1s[k])
			inisep12.append(inisep1s[k])
			minsep22.append(minsep2s[k])
			inisep22.append(inisep2s[k])
			masse2.append(masses[k])
			radialosc2.append(radialoscs[k])
			longitosc12.append(longitosc1s[k])
			longitosc22.append(longitosc2s[k])
			longitosc32.append(longitosc3s[k])

	return inisep12, minsep12, inisep22, minsep22, masse2, radialosc2, longitosc12, longitosc22, longitosc32

def categorize_by_masses_4(M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3):
	minsep1s = []
	inisep1s = []
	minsep2s = []
	inisep2s = []
	minsep3s = []
	inisep3s = []
	masses = []
	radialoscs = []
	longitoscs = []

	minsep1cur = []
	inisep1cur = []
	minsep2cur = []
	inisep2cur = []
	minsep3cur = []
	inisep3cur = []
	masscur = []
	radosccur = []
	lonosccur = []


	minsep1cur.append(min_sep1[0])
	inisep1cur.append(initial_sep1[0])
	minsep2cur.append(min_sep2[0])
	inisep2cur.append(initial_sep2[0])
	minsep3cur.append(min_sep3[0])
	inisep3cur.append(initial_sep3[0])
	masscur.append(M[0])
	radosccur.append(delta_a[0])
	lonosccur.append(delta_lambda[0])

	for i in range(1, len(min_sep1)):
		if(M[i] != M[i-1]):
			minsep1s.append(minsep1cur)
			inisep1s.append(inisep1cur)
			minsep2s.append(minsep2cur)
			inisep2s.append(inisep2cur)
			minsep3s.append(minsep3cur)
			inisep3s.append(inisep3cur)
			masses.append(masscur)
			radialoscs.append(radosccur)
			longitoscs.append(lonosccur)


			minsep1cur = []
			inisep1cur = []
			minsep2cur = []
			inisep2cur = []
			minsep3cur = []
			inisep3cur = []
			masscur = []
			radosccur = []
			lonosccur = []


		minsep1cur.append(min_sep1[i])
		inisep1cur.append(initial_sep1[i])
		minsep2cur.append(min_sep2[i])
		inisep2cur.append(initial_sep2[i])
		minsep3cur.append(min_sep3[i])
		inisep3cur.append(initial_sep3[i])
		masscur.append(M[i])
		radosccur.append(delta_a[i])
		lonosccur.append(delta_lambda[i])
	else:
		minsep1s.append(minsep1cur)
		inisep1s.append(inisep1cur)
		minsep2s.append(minsep2cur)
		inisep2s.append(inisep2cur)
		minsep3s.append(minsep3cur)
		inisep3s.append(inisep3cur)
		masses.append(masscur)
		radialoscs.append(radosccur)
		longitoscs.append(lonosccur)

	return inisep1s, minsep1s, inisep2s, minsep2s, inisep3s, minsep3s, masses, radialoscs, longitoscs 

def categorize_by_inisep_4(M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3):
	minsep1s = []
	inisep1s = []
	minsep2s = []
	inisep2s = []
	minsep3s = []
	inisep3s = []
	masses = []
	radialoscs = []
	longitoscs = []

	for j in np.linspace(28.998, 98.998 + 2.5, 30): 
		minsep1cur = []
		inisep1cur = []
		minsep2cur = []
		inisep2cur = []
		minsep3cur = []
		inisep3cur = []
		masscur = []
		radosccur = []
		lonosccur = []
		#print(j, np.round(j, 3))
		for i in range(len(min_sep1)):
			#if(abs(j - initial_sep1[i]) < 1.0):
				#print(np.round(j, 3), np.round(initial_sep1[i], 3))
			if(np.round(j, 3) == np.round(initial_sep1[i], 3)):
				minsep1cur.append(min_sep1[i])
				inisep1cur.append(str(np.round(initial_sep1[i], 3)))
				minsep2cur.append(min_sep2[i])
				inisep2cur.append(str(np.round(initial_sep2[i], 3)))
				minsep3cur.append(min_sep3[i])
				inisep3cur.append(str(np.round(initial_sep3[i], 3)))
				masscur.append(M[i])
				radosccur.append(delta_a[i])
				lonosccur.append(delta_lambda[i])
		minsep1s.append(minsep1cur)
		inisep1s.append(inisep1cur)
		minsep2s.append(minsep2cur)
		inisep2s.append(inisep2cur)
		minsep3s.append(minsep3cur)
		inisep3s.append(inisep3cur)
		masses.append(masscur)
		radialoscs.append(radosccur)
		longitoscs.append(lonosccur)

	minsep12 = []
	inisep12 = []
	minsep22 = []
	inisep22 = []
	minsep32 = []
	inisep32 = []
	masse2 = []
	radialosc2 = []
	longitosc2 = []

	for k in range(len(inisep1s)):
		if(len(inisep1s[k]) > 0):
			minsep12.append(minsep1s[k])
			inisep12.append(inisep1s[k])
			minsep22.append(minsep2s[k])
			inisep22.append(inisep2s[k])
			minsep32.append(minsep3s[k])
			inisep32.append(inisep3s[k])
			masse2.append(masses[k])
			radialosc2.append(radialoscs[k])
			longitosc2.append(longitoscs[k])

	return inisep12, minsep12, inisep22, minsep22, inisep32, minsep32, masse2, radialosc2, longitosc2	

def categorize_by_masses_5(M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3, min_sep4, initial_sep4):
	minsep1s = []
	inisep1s = []
	minsep2s = []
	inisep2s = []
	minsep3s = []
	inisep3s = []
	minsep4s = []
	inisep4s = []
	masses = []
	radialoscs = []
	longitoscs = []

	minsep1cur = []
	inisep1cur = []
	minsep2cur = []
	inisep2cur = []
	minsep3cur = []
	inisep3cur = []
	minsep4cur = []
	inisep4cur = []
	masscur = []
	radosccur = []
	lonosccur = []


	minsep1cur.append(min_sep1[0])
	inisep1cur.append(initial_sep1[0])
	minsep2cur.append(min_sep2[0])
	inisep2cur.append(initial_sep2[0])
	minsep3cur.append(min_sep3[0])
	inisep3cur.append(initial_sep3[0])
	minsep4cur.append(min_sep4[0])
	inisep4cur.append(initial_sep4[0])
	masscur.append(M[0])
	radosccur.append(delta_a[0])
	lonosccur.append(delta_lambda[0])

	for i in range(1, len(min_sep1)):
		if(M[i] != M[i-1]):
			minsep1s.append(minsep1cur)
			inisep1s.append(inisep1cur)
			minsep2s.append(minsep2cur)
			inisep2s.append(inisep2cur)
			minsep3s.append(minsep3cur)
			inisep3s.append(inisep3cur)
			minsep4s.append(minsep4cur)
			inisep4s.append(inisep4cur)
			masses.append(masscur)
			radialoscs.append(radosccur)
			longitoscs.append(lonosccur)


			minsep1cur = []
			inisep1cur = []
			minsep2cur = []
			inisep2cur = []
			minsep3cur = []
			inisep3cur = []
			minsep4cur = []
			inisep4cur = []
			masscur = []
			radosccur = []
			lonosccur = []


		minsep1cur.append(min_sep1[i])
		inisep1cur.append(initial_sep1[i])
		minsep2cur.append(min_sep2[i])
		inisep2cur.append(initial_sep2[i])
		minsep3cur.append(min_sep3[i])
		inisep3cur.append(initial_sep3[i])
		minsep4cur.append(min_sep4[i])
		inisep4cur.append(initial_sep4[i])
		masscur.append(M[i])
		radosccur.append(delta_a[i])
		lonosccur.append(delta_lambda[i])
	else:
		minsep1s.append(minsep1cur)
		inisep1s.append(inisep1cur)
		minsep2s.append(minsep2cur)
		inisep2s.append(inisep2cur)
		minsep3s.append(minsep3cur)
		inisep3s.append(inisep3cur)
		minsep4s.append(minsep4cur)
		inisep4s.append(inisep4cur)
		masses.append(masscur)
		radialoscs.append(radosccur)
		longitoscs.append(lonosccur)

	return inisep1s, minsep1s, inisep2s, minsep2s, inisep3s, minsep3s, inisep4s, minsep4s, masses, radialoscs, longitoscs 

def categorize_by_inisep_5(M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3, min_sep4, initial_sep4):
	minsep1s = []
	inisep1s = []
	minsep2s = []
	inisep2s = []
	minsep3s = []
	inisep3s = []
	minsep4s = []
	inisep4s = []
	masses = []
	radialoscs = []
	longitoscs = []
	for j in np.linspace(13.2015, 90.7015 + 2.5, 33): 
		minsep1cur = []
		inisep1cur = []
		minsep2cur = []
		inisep2cur = []
		minsep3cur = []
		inisep3cur = []
		minsep4cur = []
		inisep4cur = []
		masscur = []
		radosccur = []
		lonosccur = []
		#print(j, np.round(j, 3))
		for i in range(len(min_sep1)):
			#if(abs(j - initial_sep1[i]) < 1.0):
				#print(np.round(j, 3), np.round(initial_sep1[i], 3))
			if(np.round(j, 4) == np.round(initial_sep1[i], 4)):
				minsep1cur.append(min_sep1[i])
				inisep1cur.append(str(np.round(initial_sep1[i], 4)))
				minsep2cur.append(min_sep2[i])
				inisep2cur.append(str(np.round(initial_sep2[i], 4)))
				minsep3cur.append(min_sep3[i])
				inisep3cur.append(str(np.round(initial_sep3[i], 4)))
				minsep4cur.append(min_sep4[i])
				inisep4cur.append(str(np.round(initial_sep4[i], 4)))
				masscur.append(M[i])
				radosccur.append(delta_a[i])
				lonosccur.append(delta_lambda[i])
		minsep1s.append(minsep1cur)
		inisep1s.append(inisep1cur)
		minsep2s.append(minsep2cur)
		inisep2s.append(inisep2cur)
		minsep3s.append(minsep3cur)
		inisep3s.append(inisep3cur)
		minsep4s.append(minsep4cur)
		inisep4s.append(inisep4cur)
		masses.append(masscur)
		radialoscs.append(radosccur)
		longitoscs.append(lonosccur)

	minsep12 = []
	inisep12 = []
	minsep22 = []
	inisep22 = []
	minsep32 = []
	inisep32 = []
	minsep42 = []
	inisep42 = []
	masse2 = []
	radialosc2 = []
	longitosc2 = []

	for k in range(len(inisep1s)):
		if(len(inisep1s[k]) > 0):
			minsep12.append(minsep1s[k])
			inisep12.append(inisep1s[k])
			minsep22.append(minsep2s[k])
			inisep22.append(inisep2s[k])
			minsep32.append(minsep3s[k])
			inisep32.append(inisep3s[k])
			minsep42.append(minsep4s[k])
			inisep42.append(inisep4s[k])
			masse2.append(masses[k])
			radialosc2.append(radialoscs[k])
			longitosc2.append(longitoscs[k])

	return inisep12, minsep12, inisep22, minsep22, inisep32, minsep32, inisep42, minsep42, masse2, radialosc2, longitosc2	

# =================================================================================================================

def main(N=2):
	if(N == 2):
		M, delta_a, delta_lambda, min_sep, initial_sep = import_data_2()
		iniseps, minseps, masses, radialoscs, longitoscs = categorize_by_masses_2(M, delta_a, delta_lambda, min_sep, initial_sep)
		inisep2, minsep2, masse2, radialosc2, longitosc2 = categorize_by_inisep_2(M, delta_a, delta_lambda, min_sep, initial_sep)
		plot_delta_a_by_initial_separation_2(iniseps, radialoscs, masses)
		plot_delta_lambda_by_initial_separation_2(iniseps, longitoscs, masses)
		plot_minimum_separation_by_initial_separation_2(iniseps, minseps, masses)

		plot_delta_a_by_mass_2(inisep2, radialosc2, masse2)
		plot_delta_lambda_by_mass_2(inisep2, longitosc2, masse2)
		plot_minimum_separation_by_mass_2(inisep2, minsep2, masse2)
	elif(N == 3):
		M, delta_a, min_sep1, initial_sep1, min_sep2, initial_sep2, dlmda1, dlmda2, dlmda3 = import_data_3()
		inisep1s, minsep1s, inisep2s, minsep2s, masses, radialoscs, longitosc1s, longitosc2s, longitosc3s = categorize_by_masses_3(M, delta_a, dlmda1, dlmda2, dlmda3, min_sep1, initial_sep1, min_sep2, initial_sep2)
		inisep12, minsep12, inisep22, minsep22, masse2, radialosc2, longitosc12, longitosc22, longitosc32 = categorize_by_inisep_3(M, delta_a, dlmda1, dlmda2, dlmda3, min_sep1, initial_sep1, min_sep2, initial_sep2)
		plot_delta_a_by_initial_separation_3(inisep1s, inisep2s, radialoscs, masses)
		plot_delta_lambda_by_initial_separation_3(inisep1s, inisep2s, longitosc1s, longitosc2s, longitosc3s, masses)
		plot_minimum_separation_by_initial_separation_3(inisep1s, inisep2s, minsep1s, minsep2s, masses)

		plot_delta_a_by_mass_3(inisep12, inisep22, radialosc2, masse2)
		plot_delta_lambda_by_mass_3(inisep12, inisep22, longitosc12, longitosc22, longitosc32, masse2)
		plot_minimum_separation_by_mass_3(inisep12, minsep12, inisep22, minsep22, masse2)
	elif(N == 4):
		M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3 = import_data_4()
		inisep1s, minsep1s, inisep2s, minsep2s, inisep3s, minsep3s, masses, radialoscs, longitoscs = categorize_by_masses_4(M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3)
		inisep12, minsep12, inisep22, minsep22, inisep32, minsep32, masse2, radialosc2, longitosc2 = categorize_by_inisep_4(M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3)
		plot_delta_a_by_initial_separation_4(inisep1s, inisep2s, inisep3s, radialoscs, masses)
		plot_delta_lambda_by_initial_separation_4(inisep1s, inisep2s, inisep3s, longitoscs, masses)
		plot_minimum_separation_by_initial_separation_4(inisep1s, inisep2s, inisep3s, minsep1s, minsep2s, minsep3s, masses)

		plot_delta_a_by_mass_4(inisep12, inisep22, inisep32, radialosc2, masse2)
		plot_delta_lambda_by_mass_4(inisep12, inisep22, inisep32, longitosc2, masse2)
		plot_minimum_separation_by_mass_4(inisep12, minsep12, inisep22, minsep22, inisep32, minsep32, masse2)
	elif(N == 5):
		M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3, min_sep4, initial_sep4 = import_data_5()
		# functions that need to be written
		inisep1s, minsep1s, inisep2s, minsep2s, inisep3s, minsep3s, inisep4s, minsep4s, masses, radialoscs, longitoscs = categorize_by_masses_5(M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3, min_sep4, initial_sep4)
		inisep12, minsep12, inisep22, minsep22, inisep32, minsep32, inisep42, minsep42, masse2, radialosc2, longitosc2 = categorize_by_inisep_5(M, delta_a, delta_lambda, min_sep1, initial_sep1, min_sep2, initial_sep2, min_sep3, initial_sep3, min_sep4, initial_sep4)
		plot_delta_a_by_initial_separation_5(inisep1s, inisep2s, inisep3s, inisep4s, radialoscs, masses)
		plot_delta_lambda_by_initial_separation_5(inisep1s, inisep2s, inisep3s, inisep4s, longitoscs, masses)
		plot_minimum_separation_by_initial_separation_5(inisep1s, inisep2s, inisep3s, inisep4s, minsep1s, minsep2s, minsep3s, minsep4s, masses)

		plot_delta_a_by_mass_5(inisep12, inisep22, inisep32, inisep42, radialosc2, masse2)
		plot_delta_lambda_by_mass_5(inisep12, inisep22, inisep32, inisep42, longitosc2, masse2)
		plot_minimum_separation_by_mass_5(inisep12, minsep12, inisep22, minsep22, inisep32, minsep32, inisep42, minsep42, masse2)
	

main(3)
#main(4)
#main(5)