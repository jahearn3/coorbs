# Plot Longitude over Time
#
# Author: Joseph A'Hearn
# Created 02/14/2019
#
# This program plots the final positions of bodies 
#   in real space 
#   

import numpy as np 
import bodies as bd 
import data_loader as dl 
import plot_assistant as pa 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import constants_of_mercury6 as cm6
import useful as uf
import orbital_motion as om 
#import time 

def initialize_matrices(n_bodies, n_data):
	return np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data))

def arctangent(x, y):
	if(x == 0):
		if(y > 0):
			L = np.pi / 2
		else:
			L = (3 / 2) * np.pi
	elif(x < 0):
		L = np.arctan(y / x) + np.pi
	elif(np.arctan(y / x) < 0):
		L = np.arctan(y / x) + (2 * np.pi)
	else:
		L = np.arctan(y / x)
	return L

def initialize_four_elements(s): 	
	return s, 0, 0, 0  

def mean_motion(R_Sat, mu_Sat, J2, J4, J6, a, e, lmda, peri):
	return ((mu_Sat/(a**3.0))**0.5) * ((1.0 + (0.75 * ((R_Sat/a)**2.0) * J2) - (0.9375 * ((R_Sat/a)**4.0) * J4) + (1.09375 * ((R_Sat/a)**6.0) * J6)) + (- (0.28125 * ((R_Sat/a)**4.0) * J2**2.0) + (0.703125 * ((R_Sat/a)**6.0) * J2 * J4) + (0.2109375 * ((R_Sat/a)**6.0) * J2**3.0) + (3.0 * ((R_Sat/a)**2.0) * J2 * e**2.0)))

def horizontal_epicyclic(R_Sat, mu_Sat, J2, J4, J6, a):
	return np.sqrt(mu_Sat / (a**3)) * (1 - ((3 / 4) * J2 * ((R_Sat/a)**2)) + ((45 / 16) * J4 * ((R_Sat/a)**4)) - ((175 / 32) * J6 * ((R_Sat/a)**6)) + (- ((9 / 32) * (J2**2) * ((R_Sat/a)**4))) + ((135 / 64) * J2 * J4 * ((R_Sat/a)**6)) - ((27 / 128) * (J2**3) * ((R_Sat/a)**6)))

def vertical_epicyclic(R_Sat, mu_Sat, J2, J4, J6, a, e):
	return np.sqrt(mu_Sat / (a**3)) * (1 + ((9 / 4) * ((R_Sat/a)**2) * J2) - ((75 / 16) * ((R_Sat/a)**4) * J4) + ((245 / 32) * ((R_Sat/a)**6) * J6) + (- ((81 / 32) * ((R_Sat/a)**4) * (J2**2))) + ((675 / 64) * ((R_Sat/a)**6) * J2 * J4) + ((729 / 128) * ((R_Sat/a)**6) * (J2**3)) + (6 * ((R_Sat/a)**2) * J2 * (e**2)) - ((51 / 4) * ((R_Sat/a)**2) * J2 * (e**2)))

def all_frequencies(R_Sat, mu_Sat, J2, J4, J6, a, e, lmda, peri):
	n = mean_motion(R_Sat, mu_Sat, J2, J4, J6, a, e, lmda, peri)
	kappa = horizontal_epicyclic(R_Sat, mu_Sat, J2, J4, J6, a)
	nu = vertical_epicyclic(R_Sat, mu_Sat, J2, J4, J6, a, e)
	eta_sq = (mu_Sat / (a**3)) * (1 -  (2.0  * ((R_Sat/a)**2) * J2) + (9.375  * ((R_Sat/a)**4) * J4) - (21.875  * ((R_Sat/a)**6) * J6))
	chi_sq = (mu_Sat / (a**3)) * (1 +  (7.5  * ((R_Sat/a)**2) * J2) - (21.875 * ((R_Sat/a)**4) * J4) + (45.9375 * ((R_Sat/a)**6) * J6))
	return n, kappa, nu, eta_sq, chi_sq, (((2 * nu) + kappa) / 3), ((2 * nu) - kappa), (((2 * nu) + kappa) / 3) * ((2 * nu) - kappa)

def compute_geometric_elements(s, L, s_dot, L_dot, s_C, L_C, s_dot_C, L_dot_C, n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq):
	a = (s - s_C) / (1 - ((L_dot - L_dot_C - n) / (2 * n)))
	e = np.sqrt(((L_dot - L_dot_C - n) / (2 * n))**2 + ((s_dot - s_dot_C)/(a * kappa))**2)
	lmda = L - L_C - (2 * (n / kappa) * ((s_dot - s_dot_C) / (a * kappa)))
	if(e < 1.0E-13): # to prevent longitude of pericenter from messing things up when the orbit is nearly circular
		peri = 0
	else: 
		peri = lmda - arctangent(a * kappa * (1 - ((s - s_C) / a)), s_dot - s_dot_C)
	return a, e, lmda, peri

def compute_second_order(a, e, lmda, peri, n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq):
	s_C = a * (e**2) * ((3 / 2) * (eta_sq / (kappa**2)) - 1 - ((eta_sq / (2 * kappa**2)) * np.cos(2 * (lmda - peri)))) 
	L_C = (e**2) * ((3 / 4) + (eta_sq / (2 * (kappa**2)))) * (n / kappa) * np.sin(2 * (lmda - peri))
	s_dot_C = a * (e**2) * (eta_sq / kappa) * np.sin(2 * (lmda - peri))
	L_dot_C = ((e**2) * n * ((7 / 2) - (kappa**2 / (2 * (n**2))) - (3 * eta_sq / (kappa**2)) + (((3 / 2) + (eta_sq / (kappa**2))) * np.cos(2 * (lmda - peri))))) 
	return s_C, L_C, s_dot_C, L_dot_C

def compute_lmda_from_xyuv(x, y, u, v, AU=cm6.AU(), R_Sat=bd.R_Sat(), mu_Sat=bd.mu_Sat(), J2=bd.J2(), J4=bd.J4(), J6=bd.J6()):
	s = np.sqrt(x**2 + y**2) * AU
	L = arctangent(x, y)
	s_dot =  (( u * np.cos(L)) + (v * np.sin(L)))      * AU / 8.64E+04
	L_dot = (((-u * np.sin(L)) + (v * np.cos(L))) / s) * AU / 8.64E+04
	a, e, lmda, peri = initialize_four_elements(s)
	n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq = all_frequencies(R_Sat, mu_Sat, J2, J4, J6, a, e, lmda, peri)
	s_C, L_C, s_dot_C, L_dot_C = initialize_four_elements(0)
	a_prev = 0 
	iterate = True
	iteration = 0
	max_iterations = 150
	epsilon = 1.0E-04
	while(iterate == True):
		a, e, lmda, peri = compute_geometric_elements(s, L, s_dot, L_dot, s_C, L_C, s_dot_C, L_dot_C, n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq)
		if(e > 1):
			iterate = False
			break
		s_C, L_C, s_dot_C, L_dot_C = compute_second_order(a, e, lmda, peri, n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq)
		if(np.absolute(a - a_prev) < epsilon):
			iterate = False
		if(iteration + 1 == max_iterations):
			iterate = False
		a_prev = a
		n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq = all_frequencies(R_Sat, mu_Sat, J2, J4, J6, a, e, lmda, peri)
		iteration += 1
	if(e > 1):
		mnlng = np.nan
	else:
		mnlng = np.rad2deg(lmda) % 360
	return mnlng, e, a
	#return mnlng, e, s

def compute_lmda0(lmda, n_data):
	lmda0 = np.zeros(n_data)
	for i in range(len(lmda[0])):
		x_j = 0
		y_j = 0
		for j in range(len(lmda)):
			x_j += np.cos(np.deg2rad(lmda[j][i]))
			y_j += np.sin(np.deg2rad(lmda[j][i]))
		lmda0[i] = np.rad2deg(np.arctan2(y_j, x_j))
	return lmda0

def plot_lmda_vs_t(body_list, mass_list, t, lmda, lmda0, black=False, ncols=8, nrows=6, dpi=300):
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
	else:
		fig = plt.figure(figsize=(ncols,nrows))
	ax = plt.subplot2grid((nrows,ncols),(0, 0), rowspan=nrows, colspan=ncols - 2)
	fig, ax = pa.title_and_axes(fig, ax, 'Corotating longitude', 't [yr]', r'$\lambda - \lambda_0$' + ' [deg]', t[0], t[-1], -180, 180)
	ax.set_prop_cycle('color', plt.cm.nipy_spectral(np.linspace(0.1, 0.9, len(body_list))))
	#size = 0.001
	size = 30
	# dashed lines at stable solutions
	#ax.plot([t[0], t[-1]], [-58, -58], color='k', linestyle='dashed')
	#ax.plot([t[0], t[-1]], [-32, -32], color='k', linestyle='dashed')
	#ax.plot([t[0], t[-1]], [  0,   0], color='k', linestyle='dashed')
	#ax.plot([t[0], t[-1]], [ 29,  29], color='k', linestyle='dashed')
	#ax.plot([t[0], t[-1]], [ 62,  62], color='k', linestyle='dashed')

	ax.scatter(t, (((lmda[0] - lmda0) + 180) % 360) - 180, color='m', label=body_list[0] + '\nM = ' + str("{:.1e}".format(mass_list[0])) + ' kg', s=size)
	if(len(body_list) > 1):
		ax.scatter(t, (((lmda[1] - lmda0) + 180) % 360) - 180, color='c', label=body_list[1] + '\nM = ' + str("{:.1e}".format(mass_list[1])) + ' kg', s=size)
		ax.scatter(t, (((lmda[2] - lmda0) + 180) % 360) - 180, color='g', label=body_list[2] + '\nM = ' + str("{:.1e}".format(mass_list[2])) + ' kg', s=size)
		ax.scatter(t, (((lmda[3] - lmda0) + 180) % 360) - 180, color='darkorange', label=body_list[3] + '\nM = ' + str("{:.1e}".format(mass_list[3])) + ' kg', s=size)
		ax.scatter(t, (((lmda[4] - lmda0) + 180) % 360) - 180, color='r', label=body_list[4] + '\nM = ' + str("{:.1e}".format(mass_list[4])) + ' kg', s=size)

	exp0 = int(str("{:.1e}".format(mass_list[0])[-2:]))
	if(len(body_list) > 1):
		exp1 = int(str("{:.1e}".format(mass_list[1])[-2:]))
		exp2 = int(str("{:.1e}".format(mass_list[2])[-2:]))
		exp3 = int(str("{:.1e}".format(mass_list[3])[-2:]))
		exp4 = int(str("{:.1e}".format(mass_list[4])[-2:]))

	ax.text(t[-1] * 1.05, -80, 'Clump T'  + '\nM = ' + str("{:.1e}".format(mass_list[0])[0:3]) + r'$\times 10^{' + str(exp0) + '}$' + ' kg', color='m'         , fontsize=14, verticalalignment='center')
	if(len(body_list) > 1):
		ax.text(t[-1] * 1.05, -40, 'Clump M'  + '\nM = ' + str("{:.1e}".format(mass_list[1])[0:3]) + r'$\times 10^{' + str(exp1) + '}$' + ' kg', color='c'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05,   0, 'Clump L'  + '\nM = ' + str("{:.1e}".format(mass_list[2])[0:3]) + r'$\times 10^{' + str(exp2) + '}$' + ' kg', color='g'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05,  40, 'Clump LL' + '\nM = ' + str("{:.1e}".format(mass_list[3])[0:3]) + r'$\times 10^{' + str(exp3) + '}$' + ' kg', color='darkorange', fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05,  80, 'Object 5' + '\nM = ' + str("{:.1e}".format(mass_list[4])[0:3]) + r'$\times 10^{' + str(exp4) + '}$' + ' kg', color='r'         , fontsize=14, verticalalignment='center')

	#for i in range(len(body_list)):
	#	ax.scatter(t, (((lmda[i] - lmda0) + 180) % 360) - 180, label=body_list[i] + '\nM = ' + str("{:.1e}".format(mass_list[i])) + ' kg') 
	#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
	#plt.tight_layout()
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=14)
	pa.save_and_clear_plot(fig, ax, filename='lmda_vs_t_' + str(black), black=black, dpi=dpi)

def print_final_angular_separation(body_list, lmda, lmda0):
	print('Final angular separation between bodies: ')
	for i in range(1, len(body_list) - 1):
		print(body_list[i-1] + ' and ' + body_list[i] + ': ' + str(np.absolute(lmda[i][-1] - lmda[i-1][-1])))

def compute_deviation(lmda, lmda0):
	dev_per_body = 0
	for i in range(len(lmda)):
		variance = 0
		lmda_diff = (((lmda[i] - lmda0) + 180) % 360) - 180
		#print(np.amin(lmda_diff), np.amax(lmda_diff))
		if(i == 0):
			lmda_stationary = -58
		elif(i == 1):
			lmda_stationary = -32
		elif(i == 2):
			lmda_stationary = 0
		elif(i == 3):
			lmda_stationary = 29
		elif(i == 4):
			lmda_stationary = 62

		for j in range(len(lmda[0])):
			variance += (lmda_diff[j] - lmda_stationary)**2
		dev = np.sqrt(variance / (len(lmda[0])))
		dev_per_body += dev
	dev_per_body /= len(lmda)
	#print('standard deviation per body: ' + str(dev_per_body))
	print('deviation from stationary configuration: ' + str(dev_per_body))

def plot_s_vs_t(body_list, mass_list, t, s, black=False, ncols=8, nrows=6, dpi=100):
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
	else:
		fig = plt.figure(figsize=(ncols,nrows))
	ax = plt.subplot2grid((nrows,ncols),(0, 0), rowspan=nrows, colspan=ncols - 2)
	fig, ax = pa.title_and_axes(fig, ax, 'Radial distance from Saturn', 't [yr]', 'r', t[0], t[100])
	ax.set_prop_cycle('color', plt.cm.nipy_spectral(np.linspace(0.1, 0.9, len(body_list))))
	ax.plot(t, s[0], color='m', label=body_list[0] + '\nM = ' + str("{:.1e}".format(mass_list[0])) + ' kg')
	if(len(body_list) > 1):
		ax.plot(t, s[1], color='c', label=body_list[1] + '\nM = ' + str("{:.1e}".format(mass_list[1])) + ' kg')
		ax.plot(t, s[2], color='g', label=body_list[2] + '\nM = ' + str("{:.1e}".format(mass_list[2])) + ' kg')
		ax.plot(t, s[3], color='darkorange', label=body_list[3] + '\nM = ' + str("{:.1e}".format(mass_list[3])) + ' kg')
		ax.plot(t, s[4], color='r', label=body_list[4] + '\nM = ' + str("{:.1e}".format(mass_list[4])) + ' kg')

	exp0 = int(str("{:.1e}".format(mass_list[0])[-2:]))
	if(len(body_list) > 1):
		exp1 = int(str("{:.1e}".format(mass_list[1])[-2:]))
		exp2 = int(str("{:.1e}".format(mass_list[2])[-2:]))
		exp3 = int(str("{:.1e}".format(mass_list[3])[-2:]))
		exp4 = int(str("{:.1e}".format(mass_list[4])[-2:]))

	ymin, ymax = ax.get_ylim()

	ax.text(t[-1] * 1.05, ymax * 1 / 6, 'Clump T'  + '\nM = ' + str("{:.1e}".format(mass_list[0])[0:3]) + r'$\times 10^{' + str(exp0) + '}$' + ' kg', color='m'         , fontsize=14, verticalalignment='center')
	if(len(body_list) > 1):
		ax.text(t[-1] * 1.05, ymax * 2 / 6, 'Clump M'  + '\nM = ' + str("{:.1e}".format(mass_list[1])[0:3]) + r'$\times 10^{' + str(exp1) + '}$' + ' kg', color='c'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 3 / 6, 'Clump L'  + '\nM = ' + str("{:.1e}".format(mass_list[2])[0:3]) + r'$\times 10^{' + str(exp2) + '}$' + ' kg', color='g'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 4 / 6, 'Clump LL' + '\nM = ' + str("{:.1e}".format(mass_list[3])[0:3]) + r'$\times 10^{' + str(exp3) + '}$' + ' kg', color='darkorange', fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 5 / 6, 'Object 5' + '\nM = ' + str("{:.1e}".format(mass_list[4])[0:3]) + r'$\times 10^{' + str(exp4) + '}$' + ' kg', color='r'         , fontsize=14, verticalalignment='center')

	#for i in range(len(body_list)):
	#	ax.scatter(t, (((lmda[i] - lmda0) + 180) % 360) - 180, label=body_list[i] + '\nM = ' + str("{:.1e}".format(mass_list[i])) + ' kg') 
	#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	#ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
	plt.tight_layout()
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=14)
	pa.save_and_clear_plot(fig, ax, filename='s_vs_t_' + str(black), black=black, dpi=dpi)

def plot_a_vs_t(body_list, mass_list, t, a, black=False, ncols=8, nrows=6, dpi=100):
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
	else:
		fig = plt.figure(figsize=(ncols,nrows))
	ax = plt.subplot2grid((nrows,ncols),(0, 0), rowspan=nrows, colspan=ncols - 2)
	fig, ax = pa.title_and_axes(fig, ax, 'Semi-major axis', 't [yr]', 'a', t[0], t[-1])
	ax.set_prop_cycle('color', plt.cm.nipy_spectral(np.linspace(0.1, 0.9, len(body_list))))
	ax.plot(t, a[0], color='m', label=body_list[0] + '\nM = ' + str("{:.1e}".format(mass_list[0])) + ' kg')
	if(len(body_list) > 1):
		ax.plot(t, a[1], color='c', label=body_list[1] + '\nM = ' + str("{:.1e}".format(mass_list[1])) + ' kg')
		ax.plot(t, a[2], color='g', label=body_list[2] + '\nM = ' + str("{:.1e}".format(mass_list[2])) + ' kg')
		ax.plot(t, a[3], color='darkorange', label=body_list[3] + '\nM = ' + str("{:.1e}".format(mass_list[3])) + ' kg')
		ax.plot(t, a[4], color='r', label=body_list[4] + '\nM = ' + str("{:.1e}".format(mass_list[4])) + ' kg')

	exp0 = int(str("{:.1e}".format(mass_list[0])[-2:]))
	if(len(body_list) > 1):
		exp1 = int(str("{:.1e}".format(mass_list[1])[-2:]))
		exp2 = int(str("{:.1e}".format(mass_list[2])[-2:]))
		exp3 = int(str("{:.1e}".format(mass_list[3])[-2:]))
		exp4 = int(str("{:.1e}".format(mass_list[4])[-2:]))

	ymin, ymax = ax.get_ylim()

	ax.text(t[-1] * 1.05, ymax * 1 / 6, 'Clump T'  + '\nM = ' + str("{:.1e}".format(mass_list[0])[0:3]) + r'$\times 10^{' + str(exp0) + '}$' + ' kg', color='m'         , fontsize=14, verticalalignment='center')
	if(len(body_list) > 1):
		ax.text(t[-1] * 1.05, ymax * 2 / 6, 'Clump M'  + '\nM = ' + str("{:.1e}".format(mass_list[1])[0:3]) + r'$\times 10^{' + str(exp1) + '}$' + ' kg', color='c'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 3 / 6, 'Clump L'  + '\nM = ' + str("{:.1e}".format(mass_list[2])[0:3]) + r'$\times 10^{' + str(exp2) + '}$' + ' kg', color='g'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 4 / 6, 'Clump LL' + '\nM = ' + str("{:.1e}".format(mass_list[3])[0:3]) + r'$\times 10^{' + str(exp3) + '}$' + ' kg', color='darkorange', fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 5 / 6, 'Object 5' + '\nM = ' + str("{:.1e}".format(mass_list[4])[0:3]) + r'$\times 10^{' + str(exp4) + '}$' + ' kg', color='r'         , fontsize=14, verticalalignment='center')

	#for i in range(len(body_list)):
	#	ax.scatter(t, (((lmda[i] - lmda0) + 180) % 360) - 180, label=body_list[i] + '\nM = ' + str("{:.1e}".format(mass_list[i])) + ' kg') 
	#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	#ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
	plt.tight_layout()
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=14)
	pa.save_and_clear_plot(fig, ax, filename='a_vs_t_' + str(black), black=black, dpi=dpi)

def plot_e_vs_t(body_list, mass_list, t, e, black=False, ncols=8, nrows=6, dpi=100):
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
	else:
		fig = plt.figure(figsize=(ncols,nrows))
	ax = plt.subplot2grid((nrows,ncols),(0, 0), rowspan=nrows, colspan=ncols - 2)
	fig, ax = pa.title_and_axes(fig, ax, 'Eccentricity', 't [yr]', 'e', t[0], t[500], 0)
	ax.set_prop_cycle('color', plt.cm.nipy_spectral(np.linspace(0.1, 0.9, len(body_list))))
	ax.plot(t, e[0], color='m', label=body_list[0] + '\nM = ' + str("{:.1e}".format(mass_list[0])) + ' kg')
	if(len(body_list) > 1):
		ax.plot(t, e[1], color='c', label=body_list[1] + '\nM = ' + str("{:.1e}".format(mass_list[1])) + ' kg')
		ax.plot(t, e[2], color='g', label=body_list[2] + '\nM = ' + str("{:.1e}".format(mass_list[2])) + ' kg')
		ax.plot(t, e[3], color='darkorange', label=body_list[3] + '\nM = ' + str("{:.1e}".format(mass_list[3])) + ' kg')
		ax.plot(t, e[4], color='r', label=body_list[4] + '\nM = ' + str("{:.1e}".format(mass_list[4])) + ' kg')

	exp0 = int(str("{:.1e}".format(mass_list[0])[-2:]))
	if(len(body_list) > 1):
		exp1 = int(str("{:.1e}".format(mass_list[1])[-2:]))
		exp2 = int(str("{:.1e}".format(mass_list[2])[-2:]))
		exp3 = int(str("{:.1e}".format(mass_list[3])[-2:]))
		exp4 = int(str("{:.1e}".format(mass_list[4])[-2:]))

	ymin, ymax = ax.get_ylim()

	ax.text(t[-1] * 1.05, ymax * 1 / 6, 'Clump T'  + '\nM = ' + str("{:.1e}".format(mass_list[0])[0:3]) + r'$\times 10^{' + str(exp0) + '}$' + ' kg', color='m'         , fontsize=14, verticalalignment='center')
	if(len(body_list) > 1):
		ax.text(t[-1] * 1.05, ymax * 2 / 6, 'Clump M'  + '\nM = ' + str("{:.1e}".format(mass_list[1])[0:3]) + r'$\times 10^{' + str(exp1) + '}$' + ' kg', color='c'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 3 / 6, 'Clump L'  + '\nM = ' + str("{:.1e}".format(mass_list[2])[0:3]) + r'$\times 10^{' + str(exp2) + '}$' + ' kg', color='g'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 4 / 6, 'Clump LL' + '\nM = ' + str("{:.1e}".format(mass_list[3])[0:3]) + r'$\times 10^{' + str(exp3) + '}$' + ' kg', color='darkorange', fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 5 / 6, 'Object 5' + '\nM = ' + str("{:.1e}".format(mass_list[4])[0:3]) + r'$\times 10^{' + str(exp4) + '}$' + ' kg', color='r'         , fontsize=14, verticalalignment='center')

	#for i in range(len(body_list)):
	#	ax.scatter(t, (((lmda[i] - lmda0) + 180) % 360) - 180, label=body_list[i] + '\nM = ' + str("{:.1e}".format(mass_list[i])) + ' kg') 
	#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	#ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
	plt.tight_layout()
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=14)
	pa.save_and_clear_plot(fig, ax, filename='e_vs_t_' + str(black), black=black, dpi=dpi)

def plot_e_vs_t_geo(body_list, mass_list, t, e, black=False, ncols=8, nrows=6, dpi=100):
	e_geo = dl.geo2data(body_list[0])
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
	else:
		fig = plt.figure(figsize=(ncols,nrows))
	ax = plt.subplot2grid((nrows,ncols),(0, 0), rowspan=nrows, colspan=ncols - 2)
	fig, ax = pa.title_and_axes(fig, ax, 'Eccentricity', 't [yr]', 'e', t[0], t[500], 0, 0.1)
	ax.set_prop_cycle('color', plt.cm.nipy_spectral(np.linspace(0.1, 0.9, len(body_list))))
	ax.plot(t, e[0], color='m', label=body_list[0] + '\nM = ' + str("{:.1e}".format(mass_list[0])) + ' kg')
	ax.plot(t, e_geo, color='r')
	if(len(body_list) > 1):
		ax.plot(t, e[1], color='c', label=body_list[1] + '\nM = ' + str("{:.1e}".format(mass_list[1])) + ' kg')
		ax.plot(t, e[2], color='g', label=body_list[2] + '\nM = ' + str("{:.1e}".format(mass_list[2])) + ' kg')
		ax.plot(t, e[3], color='darkorange', label=body_list[3] + '\nM = ' + str("{:.1e}".format(mass_list[3])) + ' kg')
		ax.plot(t, e[4], color='r', label=body_list[4] + '\nM = ' + str("{:.1e}".format(mass_list[4])) + ' kg')

	exp0 = int(str("{:.1e}".format(mass_list[0])[-2:]))
	if(len(body_list) > 1):
		exp1 = int(str("{:.1e}".format(mass_list[1])[-2:]))
		exp2 = int(str("{:.1e}".format(mass_list[2])[-2:]))
		exp3 = int(str("{:.1e}".format(mass_list[3])[-2:]))
		exp4 = int(str("{:.1e}".format(mass_list[4])[-2:]))

	ymin, ymax = ax.get_ylim()

	ax.text(t[-1] * 1.05, ymax * 1 / 6, 'Clump T'  + '\nM = ' + str("{:.1e}".format(mass_list[0])[0:3]) + r'$\times 10^{' + str(exp0) + '}$' + ' kg', color='m'         , fontsize=14, verticalalignment='center')
	if(len(body_list) > 1):
		ax.text(t[-1] * 1.05, ymax * 2 / 6, 'Clump M'  + '\nM = ' + str("{:.1e}".format(mass_list[1])[0:3]) + r'$\times 10^{' + str(exp1) + '}$' + ' kg', color='c'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 3 / 6, 'Clump L'  + '\nM = ' + str("{:.1e}".format(mass_list[2])[0:3]) + r'$\times 10^{' + str(exp2) + '}$' + ' kg', color='g'         , fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 4 / 6, 'Clump LL' + '\nM = ' + str("{:.1e}".format(mass_list[3])[0:3]) + r'$\times 10^{' + str(exp3) + '}$' + ' kg', color='darkorange', fontsize=14, verticalalignment='center')
		ax.text(t[-1] * 1.05, ymax * 5 / 6, 'Object 5' + '\nM = ' + str("{:.1e}".format(mass_list[4])[0:3]) + r'$\times 10^{' + str(exp4) + '}$' + ' kg', color='r'         , fontsize=14, verticalalignment='center')

	#for i in range(len(body_list)):
	#	ax.scatter(t, (((lmda[i] - lmda0) + 180) % 360) - 180, label=body_list[i] + '\nM = ' + str("{:.1e}".format(mass_list[i])) + ' kg') 
	#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	#ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
	plt.tight_layout()
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=14)
	pa.save_and_clear_plot(fig, ax, filename='e_vs_t_geo_' + str(black), black=black, dpi=dpi)

def plot_corotating_frame(body_list, mass_list, lmda0, lmda, x, y, black=False, ncols=8, nrows=6, dpi=100):
	x1c, y1c = om.convert_to_corotating_frame(np.deg2rad(lmda0), x[0] * cm6.AU() / 6.0268E+07, y[0] * cm6.AU() / 6.0268E+07)
	x2c, y2c = om.convert_to_corotating_frame(np.deg2rad(lmda0), x[1] * cm6.AU() / 6.0268E+07, y[1] * cm6.AU() / 6.0268E+07)
	x3c, y3c = om.convert_to_corotating_frame(np.deg2rad(lmda0), x[2] * cm6.AU() / 6.0268E+07, y[2] * cm6.AU() / 6.0268E+07)
	x4c, y4c = om.convert_to_corotating_frame(np.deg2rad(lmda0), x[3] * cm6.AU() / 6.0268E+07, y[3] * cm6.AU() / 6.0268E+07)
	if(len(body_list) > 4):
		x5c, y5c = om.convert_to_corotating_frame(np.deg2rad(lmda0), x[4] * cm6.AU() / 6.0268E+07, y[4] * cm6.AU() / 6.0268E+07)
	
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
	else:
		fig = plt.figure(figsize=(ncols,nrows))
	ax = plt.subplot2grid((nrows,ncols),(0, 0), rowspan=nrows, colspan=ncols - 2)
	D = 1.3
	fig, ax = pa.title_and_axes(fig, ax, 'Corotating frame', 'x [' + r'$R_{Sat}$' + ']', 'y [' + r'$R_{Sat}$' + ']', -D, D, -D, D)
	plt.gca().set_aspect('equal', adjustable='box')
	#ax.set_prop_cycle('color', plt.cm.nipy_spectral(np.linspace(0.1, 0.9, len(body_list))))
	ax.scatter(x1c, y1c, color='m', label=body_list[0] + '\nM = ' + str("{:.1e}".format(mass_list[0])) + ' kg')
	ax.scatter(x2c, y2c, color='c', label=body_list[1] + '\nM = ' + str("{:.1e}".format(mass_list[1])) + ' kg')
	ax.scatter(x3c, y3c, color='g', label=body_list[2] + '\nM = ' + str("{:.1e}".format(mass_list[2])) + ' kg')
	ax.scatter(x4c, y4c, color='darkorange', label=body_list[3] + '\nM = ' + str("{:.1e}".format(mass_list[3])) + ' kg')
	if(len(body_list) > 4):
		ax.scatter(x5c, y5c, color='r', label=body_list[4] + '\nM = ' + str("{:.1e}".format(mass_list[4])) + ' kg')
		
	exp0 = int(str("{:.1e}".format(mass_list[0])[-2:]))
	if(len(body_list) > 1):
		exp1 = int(str("{:.1e}".format(mass_list[1])[-2:]))
		exp2 = int(str("{:.1e}".format(mass_list[2])[-2:]))
		exp3 = int(str("{:.1e}".format(mass_list[3])[-2:]))
		if(len(body_list) > 4):
			exp4 = int(str("{:.1e}".format(mass_list[4])[-2:]))

	ax.text(D * 1.05, -0.80, 'Clump T'  + '\nM = ' + str("{:.1e}".format(mass_list[0])[0:3]) + r'$\times 10^{' + str(exp0) + '}$' + ' kg', color='m'         , fontsize=14, verticalalignment='center')
	if(len(body_list) > 1):
		ax.text(D * 1.05, -0.40, 'Clump M'  + '\nM = ' + str("{:.1e}".format(mass_list[1])[0:3]) + r'$\times 10^{' + str(exp1) + '}$' + ' kg', color='c'         , fontsize=14, verticalalignment='center')
		ax.text(D * 1.05,   0, 'Clump L'  + '\nM = ' + str("{:.1e}".format(mass_list[2])[0:3]) + r'$\times 10^{' + str(exp2) + '}$' + ' kg', color='g'         , fontsize=14, verticalalignment='center')
		ax.text(D * 1.05,  0.40, 'Clump LL' + '\nM = ' + str("{:.1e}".format(mass_list[3])[0:3]) + r'$\times 10^{' + str(exp3) + '}$' + ' kg', color='darkorange', fontsize=14, verticalalignment='center')
		if(len(body_list) > 4):
			ax.text(D * 1.05,  0.80, 'Object 5' + '\nM = ' + str("{:.1e}".format(mass_list[4])[0:3]) + r'$\times 10^{' + str(exp4) + '}$' + ' kg', color='r'         , fontsize=14, verticalalignment='center')

	#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	#ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
	plt.tight_layout()
	#plt.xticks(fontsize=18)
	#plt.yticks(fontsize=14)
	pa.save_and_clear_plot(fig, ax, filename='corotating_frame_' + str(black), black=black, dpi=dpi)


def main(filename='big.in', system='', black=True, aei=True):
	#start = time.clock()
	body_list = bd.full_body_list(filename)
	mass_list = bd.full_mass_list(filename)
	# get the data from the aei files 
	t = dl.time_data(body_list[0])
	t /= 365.25

	if(aei == True):
		x, y, u, v, lmda, e, s = initialize_matrices(len(body_list), len(t))
		# compute lmda or import lmda
		# check if a geo file exists
		#if():
		for i in range(len(body_list)):
			x_i, y_i, u_i, v_i = dl.xyuv_data(body_list[i])
			#for j in range(18):
			for j in range(len(x_i)):
				x[i][j] = x_i[j]
				y[i][j] = y_i[j] 
				u[i][j] = u_i[j]
				v[i][j] = v_i[j]
				lmda[i][j], e[i][j], s[i][j] = compute_lmda_from_xyuv(x[i][j], y[i][j], u[i][j], v[i][j])	
	else:
		lmda = np.zeros((len(body_list), len(t)))
		for i in range(len(body_list)):
			lmda_i = dl.geo4data(body_list[i])
			for j in range(len(lmda_i)):
				lmda[i][j] = lmda_i[j]
	# middle longitude at each point in time
	lmda0 = compute_lmda0(lmda, len(t))
	compute_deviation(lmda, lmda0)
	plot_lmda_vs_t(body_list, mass_list, t, lmda, lmda0, black)
	#plot_corotating_frame(body_list, mass_list,lmda0, lmda, x, y)
	#plot_a_vs_t(body_list, mass_list, t, s, black)
	#plot_e_vs_t(body_list, mass_list, t, e, black)
	#plot_s_vs_t(body_list, mass_list, t, s, black)
	#plot_e_vs_t_geo(body_list, mass_list, t, e, black)
	#end = time.clock()
	#print(str(end - start) + ' seconds')
	#print_final_angular_separation(body_list, lmda, lmda0)

main(black=True, aei=True)
