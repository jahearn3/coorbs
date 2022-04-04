# Simulation Summary
#
# Author: Joseph A'Hearn
# Created 03/06/2020
#
# This program provides a summary of the simulation
#   

import numpy as np 
import bodies as bd 
import data_loader as dl 
import plot_assistant as pa 
import matplotlib.pyplot as plt 
import constants_of_mercury6 as cm6
import useful as uf
import orbital_motion as om 
import matplotlib.ticker as ticker
from sklearn.linear_model import LinearRegression
import pathlib

def initialize_matrices(n_bodies, n_data):
	return np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data))

def hill_radius(a, m, M=bd.M()):
	return (a * (m / (3 * M))**(1/3)) * 1.0E-03 # converting to km

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

def mean_motion(R, mu, J2, J4, J6, a, e, lmda, peri):
	return ((mu/(a**3.0))**0.5) * ((1.0 + (0.75 * ((R/a)**2.0) * J2) - (0.9375 * ((R/a)**4.0) * J4) + (1.09375 * ((R/a)**6.0) * J6)) + (- (0.28125 * ((R/a)**4.0) * J2**2.0) + (0.703125 * ((R/a)**6.0) * J2 * J4) + (0.2109375 * ((R/a)**6.0) * J2**3.0) + (3.0 * ((R/a)**2.0) * J2 * e**2.0)))

def horizontal_epicyclic(R, mu, J2, J4, J6, a):
	return np.sqrt(mu / (a**3)) * (1 - ((3 / 4) * J2 * ((R/a)**2)) + ((45 / 16) * J4 * ((R/a)**4)) - ((175 / 32) * J6 * ((R/a)**6)) + (- ((9 / 32) * (J2**2) * ((R/a)**4))) + ((135 / 64) * J2 * J4 * ((R/a)**6)) - ((27 / 128) * (J2**3) * ((R/a)**6)))

def vertical_epicyclic(R, mu, J2, J4, J6, a, e):
	return np.sqrt(mu / (a**3)) * (1 + ((9 / 4) * ((R/a)**2) * J2) - ((75 / 16) * ((R/a)**4) * J4) + ((245 / 32) * ((R/a)**6) * J6) + (- ((81 / 32) * ((R/a)**4) * (J2**2))) + ((675 / 64) * ((R/a)**6) * J2 * J4) + ((729 / 128) * ((R/a)**6) * (J2**3)) + (6 * ((R/a)**2) * J2 * (e**2)) - ((51 / 4) * ((R/a)**2) * J2 * (e**2)))

def all_frequencies(R, mu, J2, J4, J6, a, e, lmda, peri):
	n = mean_motion(R, mu, J2, J4, J6, a, e, lmda, peri)
	kappa = horizontal_epicyclic(R, mu, J2, J4, J6, a)
	nu = vertical_epicyclic(R, mu, J2, J4, J6, a, e)
	eta_sq = (mu / (a**3)) * (1 -  (2.0  * ((R/a)**2) * J2) + (9.375  * ((R/a)**4) * J4) - (21.875  * ((R/a)**6) * J6))
	chi_sq = (mu / (a**3)) * (1 +  (7.5  * ((R/a)**2) * J2) - (21.875 * ((R/a)**4) * J4) + (45.9375 * ((R/a)**6) * J6))
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

def compute_lmda_from_xyuv(x, y, u, v, AU=cm6.AU(), R=bd.R(), mu=bd.mu(), J2=bd.J2(), J4=bd.J4(), J6=bd.J6()):
	s = np.sqrt(x**2 + y**2) * AU
	L = arctangent(x, y)
	s_dot =  (( u * np.cos(L)) + (v * np.sin(L)))      * AU / 8.64E+04
	L_dot = (((-u * np.sin(L)) + (v * np.cos(L))) / s) * AU / 8.64E+04
	a, e, lmda, peri = initialize_four_elements(s)
	n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq = all_frequencies(R, mu, J2, J4, J6, a, e, lmda, peri)
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
		n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq = all_frequencies(R, mu, J2, J4, J6, a, e, lmda, peri)
		iteration += 1
	if(e > 1):
		mnlng = np.nan
	else:
		mnlng = np.rad2deg(lmda) % 360
	return mnlng, e, a, n

def compute_E_and_L_from_xyuv(x, y, u, v, masses, AU=cm6.AU(), G=cm6.G(), R=bd.R(), mu=bd.mu(), J2=bd.J2(), J4=bd.J4(), J6=bd.J6()):
	L = 0 
	E = 0 
	for i in range(len(masses)):
		# Angular momentum
		L += masses[i] * np.sqrt(((x[i] * v[i]) - (y[i] * u[i]))**2) * (AU**2) / 8.64E+04 # This equation is simplified because it assumes that all the angular momentum is in the z-axis
		# Energy
		v_i = np.sqrt((u[i]**2) + (v[i]**2)) * AU / 8.64E+04
		r_i = np.sqrt((x[i]**2) + (y[i]**2)) * AU
		#j2term = J2 * ((R / r_i)**2) * uf.P2(z[i] * AU / r_i) # J2, J4, and J6 are set to zero, so we can skip these lines for now
		#j4term = J4 * ((R / r_i)**4) * uf.P4(z[i] * AU / r_i)
		#j6term = J6 * ((R / r_i)**6) * uf.P6(z[i] * AU / r_i)
		#U = -(mu * masses[i] / r_i) * (1 - (j2term + j4term + j6term))
		U = -(mu * masses[i] / r_i)
		T = 0.5 * masses[i] * (v_i**2)
		# energy due to each co-orbital body
		for k in range(len(masses)):
			if(i != k):
				r_k = np.sqrt(((x[k] - x[i])**2) + ((y[k] - y[i])**2)) * AU
				U += -(G * masses[i] * masses[k] / r_k)
		E += T + U
	return E, L 

def compute_lmda0(lmda, n_data):
	lmda0 = np.zeros(n_data)
	for i in range(len(lmda[0])):
		x_j = 0
		y_j = 0
		for j in range(len(lmda)):
			x_j += np.cos(np.deg2rad(lmda[j][i]))
			y_j += np.sin(np.deg2rad(lmda[j][i]))
		lmda0[i] = np.rad2deg(np.arctan2(y_j, x_j)) 
	#print('mean lmda0: ' + str(np.mean(lmda0)))
	# attempt to prevent 180 degree jumps
	if(len(lmda) == 3):
		for i in range(len(lmda0)):
			count = 0
			while(np.absolute(((lmda[1][i] - lmda0[i] + 180) % 360) - 180) > 90):
				if(count > 10):
					print('lmda1i = ' + str(lmda[1][i]))
					print('lmda0i = ' + str(lmda0[i]))
					break
	#for i in range(1, len(lmda0)):
		#if((lmda0[i] - lmda0[i-1]) > 90):
				if((((lmda[1][i] - lmda0[i] + 180) % 360) - 180) > 90):
					lmda0[i] -= 180
		#elif((lmda0[i-1] - lmda0[i]) > 90):
				elif((((lmda[1][i] - lmda0[i] + 180) % 360) -180) > 90):
					lmda0[i] += 180
				count += 1
	#print('mean lmda0: ' + str(np.mean(lmda0)))
	return lmda0 

def compute_deviation(lmda, lmda0):
	dev_per_body = 0
	for i in range(len(lmda)):
		variance = 0
		lmda_diff = (((lmda[i] - lmda0) + 180) % 360) - 180
		if(i == 0):
			lmda_stationary = 30
		elif(i == 1):
			lmda_stationary = -30

		for j in range(len(lmda[0])):
			variance += (lmda_diff[j] - lmda_stationary)**2
		dev = np.sqrt(variance / (len(lmda[0])))
		dev_per_body += dev
	dev_per_body /= len(lmda)
	#print('standard deviation per body: ' + str(dev_per_body))
	print('deviation from stationary configuration: ' + str(dev_per_body))

def compute_a0(a, N):
	a0 = np.zeros(N)
	for i in range(N):
		a0[i] = np.mean(a[i])
	return np.mean(a0)

def compute_n0(a0, masses, M=bd.M(), G=cm6.G()):
	M_tot = M + np.sum(masses) 
	mu = G * M # This mu takes into account the masses of the bodies in orbit; the mu that is defined in bodies.py does not
	return np.rad2deg(np.sqrt((mu) / ((a0)**3))) * 8.64E+04 # deg/day

def compute_n0_from_n(n, n_data, n_bodies):
	#n0 = np.zeros(n_data)
	#for i in range(n_bodies):
	#	for j in range(n_data):
	#		n0[j] = np.sum(n[i][j]) / n_bodies
	n0 = np.mean(n, axis=0)
	return n0 * 180 * 8.64E+04 / np.pi # deg/day

def compute_lmda0_from_n0(n0, t):
	return ((n0 * t) % 360) - 180

def fix_drift(lmda, lmda0, t):
	drift = np.zeros(len(lmda))
	for j in range(len(lmda)):
		drift[j] = ((lmda[j][-1] - lmda0[-1]) - (lmda[j][0] - lmda0[0]))
	total_drift = np.mean(drift)
	#print('total_drift = ' + str(total_drift))
	for i in range(len(lmda0)):
		lmda0[i] += (total_drift * (t[i] / t[-1]))

	# unwind
	unwind = np.zeros(len(lmda))
	for j in range(len(lmda)):
		slope = np.zeros(len(t))
		for i in range(1, len(t)):
			slope[i] = (((((lmda[j][i] - lmda0[i]) + 180) % 360) - 180) - ((((lmda[j][i-1] - lmda0[i-1]) + 180) % 360) - 180)) / (t[i] - t[i-1])
			if((((((lmda[j][i] - lmda0[i]) + 180) % 360) - 180) < -170) and (((((lmda[j][i-1] - lmda0[i-1]) + 180) % 360) - 180) >  170)):
				slope[i] = slope[i-1]
			if((((((lmda[j][i] - lmda0[i]) + 180) % 360) - 180) >  170) and (((((lmda[j][i-1] - lmda0[i-1]) + 180) % 360) - 180) < -170)):
				slope[i] = slope[i-1]
		slope[0] = np.mean(slope)

		unwind[j] = np.mean(slope) * t[-1]

	drift2 = np.mean(unwind)
	for i in range(len(lmda0)):
		lmda0[i] += (drift2 * (t[i] / t[-1]))




	#if((len(lmda) % 2) != 0):
	#	if(len(lmda) == 3):
	#		j = 1
	#	elif(len(lmda) == 5):
	#		j = 2
	#	total_drift = ((((lmda[j][-1] - lmda0[-1]) + 180) % 360) - 180) - ((((lmda[j][0] - lmda0[0]) + 180) % 360) - 180) # this one works for -130
	#	if((180 < abs(total_drift) < 360) or (540 < abs(total_drift) < 720)):
	#		total_drift = ((lmda[j][-1] - lmda0[-1]) - (lmda[j][0] - lmda0[0])) # this one works for -125
	#else:
	#	if(len(lmda) == 4):
	#		total_drift1 = ((((lmda[1][-1] - lmda0[-1]) + 180) % 360) - 180) - ((((lmda[1][0] - lmda0[0]) + 180) % 360) - 180)
	#		total_drift2 = ((((lmda[2][-1] - lmda0[-1]) + 180) % 360) - 180) - ((((lmda[2][0] - lmda0[0]) + 180) % 360) - 180)
	#		total_drift = (total_drift1 + total_drift2) / 2
	#	if(total_drift > 120):
	#		total_drift = ((((lmda[1][-1] - lmda0[-1]) + 180) % 360) - 180) - ((((lmda[1][0] - lmda0[0]) + 180) % 360) - 180)
	##print('total_drift = ' + str(total_drift))
	#for i in range(len(lmda0)):
	#	lmda0[i] += (total_drift * (t[i] / t[-1]))
	
	upgrade = False
	if(upgrade == True):
		slopes = np.zeros(len(lmda))
		for i in range(len(lmda)):
			lmdafit = lmda[i] - lmda0
			#lmdafit = (((lmda[i] - lmda0) + 180) % 360) - 180
			# linear fit
			model = np.polyfit(t, lmdafit, 1)
			#print(model)
			#regression_coefficent = model[0]
			#intercept = model[1]
			#lmda0 = (regression_coefficent * t) + intercept
			predict = np.poly1d(model)
			y = predict(t)
			slopes[i] = y[-1] - y[0]
		total_drift = np.mean(slopes)
		for i in range(len(lmda0)):
			lmda0[i] += (total_drift * (t[i] / t[-1]))

	return lmda0

def fourier_transform(lmda, lmda0, t):
	black = False
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', r'$\omega$' + ' [orbits' + r'$^{-1}$' + ']', '')
	#y_fourier = np.zeros((len(lmda), len(t)))
	n_samples = 2**11
	sample_freq = np.linspace(1.0E-2, 4, n_samples)
	sample_period = np.linspace(2**5, 2**13, n_samples)
	for j in range(len(lmda)):
		y_1_omega = []
		y_2_omega = []
		amplitude_fourier = []
		period_by_fourier = []
		phi = lmda[j] - lmda0
		for i in range(n_samples):
			y_1_omega.append(np.sum(phi * np.cos(sample_freq[i] * t)))
			y_2_omega.append(np.sum(phi * np.sin(sample_freq[i] * t)))
			#y_1_omega.append(np.sum(phi * np.cos(2 * np.pi * t / sample_period[i])))
			#y_2_omega.append(np.sum(phi * np.sin(2 * np.pi * t / sample_period[i])))
			amplitude_fourier.append(2 * np.sqrt(y_1_omega[i]**2 + y_2_omega[i]**2) / len(t))
			#period_by_fourier.append(2 * np.pi / sample_freq[i])
			#y_fourier[j][i] = 2 * np.sqrt(y_1_omega[i]**2 + y_2_omega[i]**2) / len(t)
		#print(np.amax(amplitude_fourier))
		#print(sample_period[np.argmax(amplitude_fourier)])
		#libration_period = np.amax(freq)
		#print(libration_period)
		ax.plot(sample_freq, amplitude_fourier, color=assigned_color(j))
		p = find_peaks(amplitude_fourier, 0.5 * np.amax(amplitude_fourier))
		ax.set_xscale('log')
		for i in range(len(p)):
			ax.plot([sample_freq[p[i]], sample_freq[p[i]]], [0, amplitude_fourier[p[i]] / np.amax(amplitude_fourier)], linestyle='dashed', label=str(np.around(2 * np.pi / sample_freq[p[i]], 1)) + ' orbits')
		#ax.plot(sample_period, amplitude_fourier, color=assigned_color(j))
	ax.legend(loc=0, fontsize=14)
	pa.save_and_clear_plot(fig, ax, filename='fourier_' + str(len(lmda)), black=black)

def find_peaks(x, h):
	p = []
	for i in range(1, len(x) - 1):
		if(x[i] > h):
			if(x[i] > x[i-1]) and (x[i] > x[i+1]):
				p.append(i)
	return p

def compute_delta_a(a, a0, R_H, units='R_H'):
	max_a = np.zeros(len(a))
	min_a = np.zeros(len(a))
	for i in range(len(a)):
		max_a[i] = np.amax(a[i])
		min_a[i] = np.amin(a[i])
	delta_a = np.amax(np.fmax(max_a - a0, a0 - min_a)) * 1.0E-03
	if(units == 'R_H'):
		return delta_a / R_H # units of Hill radii
	else:
		return delta_a # units of km

def compute_delta_lambda(lmda, lmda0):
	max_l = np.zeros(len(lmda))
	min_l = np.zeros(len(lmda))
	for i in range(len(lmda)):
		maxis = np.amax((((lmda[i] - lmda0) + 180) % 360) - 180)
		minis = np.amin((((lmda[i] - lmda0) + 180) % 360) - 180)
		if(maxis < 0):
			maxis += 360
		elif(maxis > 360):
			maxis -= 360
		if(maxis > 180):
			maxis = 180 - maxis
		if(minis < 0):
			minis += 360
		elif(minis > 360):
			minis -= 360
		if(minis > 180):
			minis = 180 - minis
		max_l[i] = maxis
		min_l[i] = minis
	delta_lambda = np.sqrt(np.amax(np.fmax(max_l**2, min_l**2)))
	if(delta_lambda < 0):
		delta_lambda += 360
	elif(delta_lambda > 360):
		delta_lambda = 360
	return delta_lambda

def compute_delta_e(e):
	max_e = np.zeros(len(e))
	for i in range(len(e)):
		max_e[i] = np.amax(e[i])
	return np.amax(max_e)

def compute_minimum_separation_deg(lmda, lmda0): # this function is deprecated; an improved version is computed within the plot functions
	#mins = np.zeros(len(lmda))
	#argmins = np.zeros(len(lmda))
	#for i in range(len(lmda)):
	#	for j in range(len(lmda)):
	#		if(i < j):
	#			sep = ((((lmda[i] - lmda0) + 180) % 360) - 180) - ((((lmda[j] - lmda0) + 180) % 360) - 180)
	#			mins[i] = np.amin(sep)
	#			argmins[i] = int(np.argmin(sep))
	mins = np.zeros(len(lmda) - 1)
	argmins = np.zeros(len(lmda) - 1)
	for i in range(len(lmda) - 1):
		sep = abs(((((lmda[i] - lmda0) + 180) % 360) - 180) - ((((lmda[i+1] - lmda0) + 180) % 360) - 180)) % 360
		for j in range(len(sep)):
			if(sep[j] > 180):
				sep[j] = 180 - sep[j]
		mins[i] = np.amin(sep)
		argmins[i] = int(np.argmin(sep))
	return mins, argmins
	

def compute_initial_separation_deg(lmda, lmda0):
	#initial_sep_deg = np.zeros(len(lmda))
	initial_sep_deg = np.zeros(len(lmda) - 1)
	#for i in range(len(lmda)):
	#	for j in range(len(lmda)):
	#		if(i < j):
	#			initial_sep_deg[i] = ((((lmda[i][0] - lmda0[0]) + 180) % 360) - 180) - ((((lmda[j][0] - lmda0[0]) + 180) % 360) - 180)
	for i in range(len(lmda) - 1):
		initial_sep_deg[i] = abs(((((lmda[i][0] - lmda0[0]) + 180) % 360) - 180) - ((((lmda[i+1][0] - lmda0[0]) + 180) % 360) - 180)) % 360
		if(initial_sep_deg[i] > 180):
			initial_sep_deg[i] = abs(180 - initial_sep_deg[i])
		if(initial_sep_deg[i] < 0):
			initial_sep_deg[i] = abs(initial_sep_deg[i])

	return initial_sep_deg

def plot_semimajor_axis(fig, bodies, a, a0, t, delta_a, R_H, units='R_H', t_units='orbits', nrows=1, ncols=1, initrow=0, initcol=0, rs=1, cs=1):
	ax1 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	if(units == 'R_H'):
		fig, ax1 = pa.title_and_axes(fig, ax1, '', 't [yr]', r'$a - a_0$' + ' [' + r'$R_H$' + ']', t[0], t[-1], -1.1*delta_a, 1.1*delta_a)
	else:
		fig, ax1 = pa.title_and_axes(fig, ax1, '', 't [yr]', r'$a - a_0$' + ' [km]', t[0], t[-1], -1.1*delta_a, 1.1*delta_a)
	if(t_units == 'orbits'):
		ax1.set_xlabel('t [orbits]', fontsize=20)
	for i in range(len(bodies)):
		if(units == 'R_H'):
			ax1.plot(t, (a[i] - a0) * 1.0E-03 / R_H, color=assigned_color(i), label=uf.decapitalize(bodies[i]), linewidth=3)
		else:
			ax1.plot(t, (a[i] - a0) * 1.0E-03, color=assigned_color(i), label=uf.decapitalize(bodies[i]), linewidth=3)
	start, end = ax1.get_xlim()
	ax1.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	ax1.legend(loc=4, prop={'size':24})
	plt.setp(ax1.get_xticklabels(), visible=False) # make these tick labels invisible
	ax1.set_xlabel('')
	#ax1.text(0.5, 0.5, str(np.round(a0 * 1.0E-03, 0)) + ' km', fontsize=16, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	return fig, ax1

def plot_corotating_longitude(fig, bodies, lmda, lmda0, delta_lambda, min_sep_deg, min_sep_idx, t, t_units, nrows, ncols, initrow, initcol, rs, cs):
	#if(len(lmda) == 3):
	#	j = 1
	#elif(len(lmda) == 5):
	#	j = 2
	#if(len(lmda) == 4):
		#lmdafi0 = lmda[0] - lmda0 - 180
		#lmdafi1 = lmda[1] - lmda0 - 180
		#lmdafi2 = lmda[2] - lmda0 - 180
		#lmdafi3 = lmda[4] - lmda0 - 180
		#lmdafit = (lmdafi1 + lmdafi2) / 2
		#lmdafit = (((lmda[1] - lmda0) + 180) % 360) - 180
	#total_drift = ((((lmda[j][-1] - lmda0[-1]) + 180) % 360) - 180) - ((((lmda[j][0] - lmda0[0]) + 180) % 360) - 180) # this one works for -130
	##for i in range(len(lmda0)):
	##	lmda0[i] += (lmdafit[-1] * (t[i] / t[-1]))
	#
	# linear fit
	#model = np.polyfit(t, lmdafit, 1)
	#print(model)
	#regression_coefficent = model[0]
	#intercept = model[1]
	##lmda0 = (regression_coefficent * t) + intercept
	#predict = np.poly1d(model)
	#y = predict(t)
	#total_drift = y[-1] - y[0]
	#for i in range(len(lmda0)):
	#	lmda0[i] += (total_drift * (t[i] / t[-1]))



	#t0 = t.reshape((-1, 1))
	#lrmodel = LinearRegression().fit(t0, lmda0)
	#r_sq = lrmodel.score(t0, lmda0)
	#print('coefficient of determination: ' + str(r_sq))
	#slope = lrmodel.coef_[0] 
	#print('slope: ' + str(slope))
	##y_pred = lrmodel.intercept_ + np.sum(lrmodel.coef_ * t0, axis=1)
	#y_pred = -slope * t


	ax3 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax3 = pa.title_and_axes(fig, ax3, '', 't [yr]', r'$\lambda - \lambda_0$' + ' [deg]', t[0], t[-1], -180, 180)
	if(t_units == 'orbits'):
		ax3.set_xlabel('t [orbits]', fontsize=20)
	for i in range(len(bodies)):
		#ax3.scatter(t, (((lmda[i] - lmda00) + 180) % 360) - 180, color=assigned_color(i), label=uf.decapitalize(bodies[i]))
		ax3.scatter(t, (((lmda[i] - lmda0) + 180) % 360) - 180, color=assigned_color(i), label=uf.decapitalize(bodies[i]))
	for i in range(len(bodies) - 1):
		min_sep_t = t[int(min_sep_idx[i])]
		min_sep_y1 = (((lmda[ i ][int(min_sep_idx[i])] - lmda0[int(min_sep_idx[i])]) + 180) % 360) - 180
		min_sep_y2 = (((lmda[i+1][int(min_sep_idx[i])] - lmda0[int(min_sep_idx[i])]) + 180) % 360) - 180
		ax3.vlines(x=min_sep_t, ymin=min_sep_y1, ymax=min_sep_y2, linestyles='dashed')
		ax3.text(t[-1] / 2, (min_sep_y1 + min_sep_y2) / 2, 'minimum separation: ' + str(np.round(min_sep_deg[i], 1)) + ' deg', fontsize=24, horizontalalignment='center', verticalalignment='center')
	ax3.set_yticks([-180, -120, -60, 0, 60, 120, 180])
	#ax3.text(t[0], 170, r'$n_0 = $' + str(n0) + 'what units?')
	#start, end = ax3.get_xlim()
	#ax3.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	#ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	#ax2.legend(loc=2, prop={'size':14})
	#plt.setp(ax10.get_xticklabels(), visible=False) # make these tick labels invisible
	#ax1.set_xlabel('')		
	ax3.spines['top'].set_visible(False)
	ax3.spines['right'].set_visible(False)

	#ax3.plot(t, y_pred, c='r')
	#ax3.plot(t, y, c='r')

	return fig, ax3

def plot_corotating_frame(fig, bodies, a, a0, lmda, lmda0, nrows, ncols, initrow, initcol, rs, cs, AU=cm6.AU()):
	ax2 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax2 = pa.title_and_axes(fig, ax2, '', 'x [' + r'$a_0$' + ']', 'y [' + r'$a_0$' + ']', -1.2, 1.2, -1.2, 1.2)
	plt.gca().set_aspect('equal', adjustable='box')

	# circle for orbit
	angles = np.linspace(0, 2 * np.pi, 1000)
	ax2.scatter(np.cos(angles), np.sin(angles), s=0.1, color='gray')

	maxradexc = 0
	for i in range(len(bodies)):
		if(np.amax(((a[i] - a0) / a0) > maxradexc)):
			maxradexc = np.amax((a[i] - a0) / a0)
	print('max radial exc: ' + str(maxradexc))
	#if(maxradexc < 0.0001):
	#	enhancement_factor = 500
	#elif(maxradexc < 0.006):
	#	enhancement_factor = 400
	#elif(maxradexc < 0.006):
	#	enhancement_factor = 300
	#elif(maxradexc < 0.006):
	#	enhancement_factor = 200
	#elif(maxradexc < 0.002):
	#	enhancement_factor = 100
	#elif(maxradexc < 0.006):
	#	enhancement_factor = 75
	#else:
	#	enhancement_factor = 50
	enhancement_factor = 5 * np.round(100 * maxradexc, decimals=0)
	for i in range(len(bodies)):
		theta = np.deg2rad((((lmda[i] - lmda0) + 180) % 360) - 180)
		radexc = (a[i] - a0) / a0
		
		#print('max radial excursion: ' + str(np.amax(radexc)))
		ax2.scatter(((a[i] * np.cos(theta)) / a0) + (enhancement_factor * radexc * np.cos(theta)), ((a[i] * np.sin(theta)) / a0) + (enhancement_factor * radexc * np.sin(theta)), s=10, color=assigned_color(i))
		# angles for initial positions
		ax2.plot([0,(a[i][0] * np.cos(theta[0])) / a0], [0,(a[i][0] * np.sin(theta[0])) / a0], color='gray')
	# dot at origin
	ax2.scatter(0, 0, 500, 'k', zorder=4)
	ax2.text(0, 0.5, 'radial excursions \nenhanced ' + str(enhancement_factor) + 'x', fontsize=24, horizontalalignment='center', verticalalignment='center')
	return fig, ax2

def plot_eccentricity(fig, bodies, e, delta_e, t, t_units, nrows, ncols, initrow, initcol, rs, cs):
	ax2 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax2 = pa.title_and_axes(fig, ax2, '', 't [yr]', r'$e$', t[0], t[-1], 0, 1.1*delta_e)
	if(t_units == 'orbits'):
		ax2.set_xlabel('t [orbits]', fontsize=20)
	for i in range(len(bodies)):
		ax2.scatter(t, e[i], color=assigned_color(i), label=uf.decapitalize(bodies[i]), s=1)
	start, end = ax2.get_xlim()
	ax2.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	#ax3.legend(loc=0, prop={'size':14})
	plt.setp(ax2.get_xticklabels(), visible=False) # make these tick labels invisible
	ax2.set_xlabel('')
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	return fig, ax2

def plot_phase_space(fig, bodies, a, a0, delta_a, lmda, lmda0, delta_lambda, R_H, units, nrows, ncols, initrow, initcol, rs, cs):
	ax4 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax4 = pa.title_and_axes(fig, ax4, '', r'$a - a_0$' + ' [km]', r'$\lambda - \lambda_0$ + 180 % 360 - 180' + ' [deg]', xmin=-1.1*delta_a, xmax=1.1*delta_a, ymin=-1.1*delta_lambda, ymax=1.1*delta_lambda)
	# error bars to account for eccentricity
	#xerri = a[i][int((begidx + endidx) / 2)] * e[i][int((begidx + endidx) / 2)] * 1.0E-03
	#xerrj = a[j][int((begidx + endidx) / 2)] * e[j][int((begidx + endidx) / 2)] * 1.0E-03
	#print('xerri: ' + str(xerri))
	#print('xerrj: ' + str(xerrj))
	#
	#ax2.errorbar((a[j][int((begidx + endidx) / 2)] - a_cer[int((begidx + endidx) / 2)]) * 1.0E-03, ((lmda[j][int((begidx + endidx) / 2)] - lambda_cer[int((begidx + endidx) / 2)] + 180) % 360) - 180, xerr=xerrj, color='c', capsize=10, capthick=2)
	#ax2.errorbar((a[i][int((begidx + endidx) / 2)] - a_cer[int((begidx + endidx) / 2)]) * 1.0E-03, ((lmda[i][int((begidx + endidx) / 2)] - lambda_cer[int((begidx + endidx) / 2)] + 180) % 360) - 180, xerr=xerri, color='r', capsize=10, capthick=2)

	for i in range(len(bodies)):
		ax4.scatter((a[i] - a0) * 1.0E-03, ((lmda[i] - lmda0 + 180) % 360) - 180, label=uf.decapitalize(bodies[i]), s=180)
	
	#ax2.legend(loc=0, prop={'size':14})
	ax4.spines['top'].set_visible(False)
	ax4.spines['right'].set_visible(False)
	return fig, ax4

def assigned_color(i):
	if(i == 0):
		return 'b'
	if(i == 1):
		return 'orange'
	if(i == 2):
		return 'g'
	if(i == 3):
		return 'm'
	if(i == 4):
		return 'r'

def assigned_color2(i):
	if(i == 0):
		return 'c'
	if(i == 1):
		return 'mediumpurple'
	if(i == 2):
		return 'olive'
	if(i == 3):
		return 'saddlebrown'
	if(i == 4):
		return 'k'

def determine_xaxis(units, x, R_H):
	if(units == 'R_H'):
		xlabel = r'$a - a_0$' + ' [' + r'$R_H$' + ']'
		x /= R_H
	else:
		xlabel = r'$a - a_0$' + ' [km]'
	return x, xlabel

def determine_yaxis(lmdai, lmda0, margin=0.1):
	# determining y limits and offsets
	maxy = np.amax(((lmdai - lmda0 + 180) % 360) - 180)
	miny = np.amin(((lmdai - lmda0 + 180) % 360) - 180)
	diff = maxy - miny
	maxy360 = np.amax(((lmdai - lmda0) % 360))
	miny360 = np.amin(((lmdai - lmda0) % 360))
	diff360 = maxy360 - miny360
	if(diff <= diff360):
		ymin = miny - (diff * margin)
		if(ymin < -180):
			ymin = -180
		ymax = maxy + (diff * margin)
		if(ymax > 180):
			ymax = 180
		y = ((lmdai - lmda0 + 180) % 360) - 180
	else: 
		ymin = miny360 - (diff * margin)
		#if(ymin < 7):
		#	ymin = 0
		ymax = maxy360 + (diff * margin)
		if(ymax > 360):
			ymax = 360
		diff = diff360
		y = ((lmdai - lmda0) % 360)
	#print('ymin: ' + str(ymin))	
	#print('ymax: ' + str(ymax))
	#print('diff: ' + str(diff))
	return y, ymin, ymax, diff

def plot_phase_space_separated(fig, bodies, i, a, a0, delta_a, lmda, lmda0, delta_lambda, R_H, units, nrows, ncols, initrow, initcol, rs, cs):
	ax4 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	x, xlabel = determine_xaxis(units, (a[i] - a0) * 1.0E-03, R_H)
	ylabel = r'$\lambda - \lambda_0$' + ' [deg]'
	y, ymin, ymax, diff = determine_yaxis(lmda[i], lmda0)
	fig, ax4 = pa.title_and_axes(fig, ax4, '', xlabel, ylabel, xmin=-1.1*delta_a, xmax=1.1*delta_a, ymin=ymin, ymax=ymax)
	ax4.scatter(x, y, color=assigned_color(i), label=uf.decapitalize(bodies[i]), s=180)
	
	if((i + 1) != len(lmda)):
		plt.setp(ax4.get_xticklabels(), visible=False) # make these tick labels invisible
		ax4.set_xlabel('')
	start, end = ax4.get_ylim()
	
	# crosshairs at equilibrium points
	s = 500
	c = 'k'
	style = '+'
	if((len(lmda) % 2) != 0):
		if(start < 0 < end):
			ax4.scatter(0, 0,       s, c, marker=style) 
	if(len(lmda) == 3): 
		if(start < -47.361 < end):	
			ax4.scatter(0, -47.361, s, c, marker=style)
		if(start <  47.361 < end):	
			ax4.scatter(0,  47.361, s, c, marker=style)
		if(start < 312.639 < end):
			ax4.scatter(0, 312.639, s, c, marker=style)
	elif(len(lmda) == 4): 
		if(start < -18.678 < end):
			ax4.scatter(0, -18.678, s, c, marker=style)
		if(start <  18.678 < end):
			ax4.scatter(0,  18.678, s, c, marker=style)
		if(start < -60.176 < end):
			ax4.scatter(0, -60.176, s, c, marker=style)
		if(start <  60.176 < end):
			ax4.scatter(0,  60.176, s, c, marker=style)
		if(start < 341.322 < end):
			ax4.scatter(0, 341.322, s, c, marker=style)
		if(start < 299.824 < end):
			ax4.scatter(0, 299.824, s, c, marker=style)
	elif(len(lmda) == 5): 
		if(start < -32.6600 < end):
			ax4.scatter(0, -32.6600, s, c, marker=style)
		if(start <  32.6600 < end):
			ax4.scatter(0,  32.6600, s, c, marker=style)
		if(start < -70.8615 < end):
			ax4.scatter(0, -70.8615, s, c, marker=style)
		if(start <  70.8615 < end):
			ax4.scatter(0,  70.8615, s, c, marker=style)
		if(start < 327.3400 < end):
			ax4.scatter(0, 327.3400, s, c, marker=style)
		if(start < 289.1385 < end):
			ax4.scatter(0, 289.1385, s, c, marker=style)
	#ax2.legend(loc=0, prop={'size':14})
	if(diff > 120):
		ax4.yaxis.set_major_locator(ticker.MultipleLocator(60))
	elif(diff > 60):
		ax4.yaxis.set_major_locator(ticker.MultipleLocator(30))
	elif(diff > 30):
		ax4.yaxis.set_major_locator(ticker.MultipleLocator(15))
	elif(diff > 20):
		ax4.yaxis.set_major_locator(ticker.MultipleLocator(10))
	elif(diff > 10):
		ax4.yaxis.set_major_locator(ticker.MultipleLocator(5))
	elif(diff > 5):
		ax4.yaxis.set_major_locator(ticker.MultipleLocator(2))
	elif(diff > 2): 
		ax4.yaxis.set_major_locator(ticker.MultipleLocator(1))
	ax4.spines['top'].set_visible(False)
	ax4.spines['right'].set_visible(False)
	
	return fig, ax4, diff

def plot_longitudinal_differences_separated(fig, bodies, i, lmda, t, t_units, nrows, ncols, initrow, initcol, rs, cs):
	ax5 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	
	if(i == (len(lmda) - 1)):
		ylabelprefix = r'$\lambda_' + str(i+1) + ' - ' + '\lambda_1$'
		j = 0
	else:
		ylabelprefix = r'$\lambda_' + str(i+1) + ' - ' + '\lambda_' + str(i+2) + '$'
		j = i + 1
	ylabel = ylabelprefix + ' [deg]'
	y, ymin, ymax, diff = determine_yaxis(lmda[i], lmda[j])
	if((i + 1) != len(lmda)):
		plt.setp(ax5.get_xticklabels(), visible=False) # make these tick labels invisible
		xlabel = ''
	elif(t_units == 'orbits'):
		xlabel = 't [orbits]'
	else:
		xlabel = 't [yr]'	
	fig, ax5 = pa.title_and_axes(fig, ax5, '', xlabel, ylabel, t[0], t[-1], ymin, ymax)	
	ax5.scatter(t, y, color=assigned_color2(i), label=uf.decapitalize(bodies[i]), s=180)
	min_sep =   np.amin(np.sqrt((((lmda[i] - lmda[j] + 180) % 360) - 180)**2))
	idx = int(np.argmin(np.sqrt((((lmda[i] - lmda[j] + 180) % 360) - 180)**2)))
	
	if(diff > 120):
		ax5.yaxis.set_major_locator(ticker.MultipleLocator(60))
	elif(diff > 60):
		ax5.yaxis.set_major_locator(ticker.MultipleLocator(30))
	elif(diff > 30):
		ax5.yaxis.set_major_locator(ticker.MultipleLocator(15))
	elif(diff > 20):
		ax5.yaxis.set_major_locator(ticker.MultipleLocator(10))
	elif(diff > 10):
		ax5.yaxis.set_major_locator(ticker.MultipleLocator(5))
	elif(diff > 5):
		ax5.yaxis.set_major_locator(ticker.MultipleLocator(2))
	elif(diff > 2): 
		ax5.yaxis.set_major_locator(ticker.MultipleLocator(1))
	
	#ax5.set_yticks([-180, -120, -60, 0, 60, 120, 180])
	#ax2.legend(loc=0, prop={'size':14})
	ax5.spines['top'].set_visible(False)
	ax5.spines['right'].set_visible(False)
	#start, end = ax5.get_ylim()
	#diff = end - start 
	return fig, ax5, min_sep, idx, ymax - ymin 

def plot_longitude_histograms(lmda, input_from_file=True, black=False): # so far this is for a single simulation
	span = 360 # or 360
	n_bins = int(span / 5)
	if(input_from_file == False):
		longitudes = []
		for i in range(len(lmda)):
			for j in range(len(lmda)):
				if(i < j):
					y, ymin, ymax, diff = determine_yaxis(lmda[i], lmda[j])
					for k in range(len(y)):
						longitudes.append(y[k])
						with open('longitudes.txt', 'a') as myfile:
							myfile.write(str(y[k]))
							myfile.write("\n")
	else:
		# lmda needs to already contain all the longitudes
		longitudes = lmda 
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(8,8))
	else:
		fig = plt.figure(figsize=(8,8))
	ax = plt.subplot2grid((8,8),(0, 0), rowspan=8, colspan=8)
	fig, ax = pa.title_and_axes(fig, ax, 'Longitudinal separations', r'$\lambda_i - \lambda_j$' + ' [deg]', 'Sample Fraction')
	#fig, ax = pa.title_and_axes(fig, ax, 'Longitudinal separations', r'$\lambda_i - \lambda_j$' + ' [deg]', 'Sample Fraction', 0, 360, 0)
	#ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
	ax.hist(longitudes, bins=n_bins, weights=np.ones(len(longitudes)) / len(longitudes), color='g')
	#print(len(longitudes))
	#print(np.amin(longitudes), np.amax(longitudes))
	#plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='longitude_histograms', black=black)
	
def longitude_histogram_from_file(filename='longitudes.txt'):
	data = np.loadtxt(filename)
	lmda = data[:]
	a_plot(['bodyA', 'bodyB'], [2.2, 2.1], 2, [0], 0.2, 1)
	plot_longitude_histograms(lmda)


def summary_plot(bodies, t, a, a0, delta_a, R_H, e, delta_e, lmda, lmda0, delta_lambda, x, y, black=False, units='R_H', t_units='orbits', nrows=60, ncols=27):
	margin = 2
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(ncols,ncols))
	else:
		fig = plt.figure(figsize=(ncols,ncols))
	fig, ax1 = plot_semimajor_axis(fig, bodies, a, a0, t, delta_a, R_H, units, t_units, nrows, ncols, 0, 0, 20 - (2 * margin), 9 - margin)
	#fig, ax2 = plot_eccentricity(fig, bodies, e, delta_e, t, t_units, nrows, ncols, 20, 0, 20, 9)
	#fig, ax2 = plot_corotating_frame(fig, bodies, x, y, a0, lmda, lmda0, nrows, ncols, 20, 0, 20, 9)
	fig, ax2 = plot_corotating_frame(fig, bodies, a, a0, lmda, lmda0, nrows, ncols, 40, 0, 20 - (2 * margin), 9 - margin)
	N = len(lmda)
	delta_lambda_max = np.zeros(N)
	min_sep_deg = np.zeros(N)
	min_sep_idx = np.zeros(N)
	for i in range(N):
		fig, ax4, diff = plot_phase_space_separated(fig, bodies, i, a, a0, delta_a, lmda, lmda0, delta_lambda, R_H, units, nrows, ncols, i * int(60 / N), 9, int(60 / N) - (2 * margin), 9 - margin)
		#if(diff > delta_lambda_max):
		#	delta_lambda_max = diff
		fig, ax5, min_sep, idx, dlmax = plot_longitudinal_differences_separated(fig, bodies, i, lmda, t, t_units, nrows, ncols, i * int(60 / N), 18, int(60 / N) - (2 * margin), 9 - margin)
		min_sep_deg[i] = min_sep
		min_sep_idx[i] = idx
		delta_lambda_max[i] = dlmax 
	#fig, ax3 = plot_corotating_longitude(fig, bodies, lmda, lmda0, delta_lambda, min_sep_deg, min_sep_idx, t, t_units, nrows, ncols, 40, 0, 20, 9)
	fig, ax3 = plot_corotating_longitude(fig, bodies, lmda, lmda0, delta_lambda, min_sep_deg, min_sep_idx, t, t_units, nrows, ncols, 20, 0, 20 - (2 * margin), 9 - margin)
	#plt.tight_layout()
	ax1.figure.savefig('simulation_summary.png')
	plt.clf()
	return delta_lambda_max, min_sep_deg

def record_data(masses, delta_a, lmda_amp, min_sep_deg, initial_sep_deg, E, L, units):
	if(units != 'R_H'):
		units = 'km'
	if(pathlib.Path('stats_tot.txt').exists == False):
	    with open('stats_tot.txt', 'a') as myfile:
	    	myfile.write('M (kg)\t\t\t') 
	    	myfile.write('delta a (' + units + ')\t')
	    	#myfile.write('delta lambda (deg)\t')
	    	for i in range(len(masses) - 1):
	    		myfile.write('min sep '     + str(i+1) + ' (deg)\t\t')
	    		myfile.write('initial sep ' + str(i+1) + ' (deg)\t\t')
	    	myfile.write('E (J)\t')
	    	myfile.write('L (kg m^2/s))\t')
	    	myfile.write('\n')
	with open('stats_tot.txt', 'a') as myfile:
		myfile.write("%.2e" % masses[0])
		myfile.write("\t\t%.4f" % delta_a)
		for i in range(len(masses) - 1):
			myfile.write("\t\t\t\t%.4f" % min_sep_deg[i])
			myfile.write("\t\t\t\t%.4f" % initial_sep_deg[i])
		myfile.write("\t\t%.4f" % E)
		myfile.write("\t\t%.4f" % L)
		#for i in range(len(masses)):
		#	myfile.write("\t\t\t%.4f" % lmda_amp[i])
		myfile.write("\n")

def a_plot(bodies, a, a0, t, delta_a, R_H, units='R_H', t_units='orbits', black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = plot_semimajor_axis(fig, bodies, a, a0, t, delta_a, R_H, units, t_units)
	pa.save_and_clear_plot(fig, ax, filename='a_vs_t', black=black) 

def E_plot(bodies, E, t, t_units='orbits', black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, 'Energy', 'E', 't (' + t_units + ')', t[0], t[-1])
	ax.plot(t, E)
	pa.save_and_clear_plot(fig, ax, filename='E_vs_t', black=black) 

def main(filename='big.in', black=True, units='R_H', t_units='orbits'):
	bodies = bd.full_body_list(filename)
	masses = bd.full_mass_list(filename)
	# get the data from the aei files 
	t = dl.time_data(bodies[0])
	x, y, u, v, lmda, e, a, n = initialize_matrices(len(bodies), len(t))

	for i in range(len(bodies)):
		x_i, y_i, u_i, v_i = dl.xyuv_data(bodies[i])
		for j in range(len(x_i)):
			x[i][j] = x_i[j]
			y[i][j] = y_i[j] 
			u[i][j] = u_i[j]
			v[i][j] = v_i[j]
			lmda[i][j], e[i][j], a[i][j], n[i][j] = compute_lmda_from_xyuv(x[i][j], y[i][j], u[i][j], v[i][j])
	E, L = compute_E_and_L_from_xyuv(x, y, u, v, masses)
	# middle longitude at each point in time
	#lmda0 = compute_lmda0(lmda, len(t))
	a0 = compute_a0(a, len(bodies))
	#n00 = compute_n0(a0, masses) # units of deg/day
	#print('n0 = ' + str(n00) + ' deg/day')
	n0 = compute_n0_from_n(n, len(t), len(bodies)) # units of deg/day
	print('mean n0 = ' + str(np.mean(n0)) + ' deg/day')
	lmda0 = compute_lmda0_from_n0(n0, t)
	lmda0 = fix_drift(lmda, lmda0, t)
	# NEED TO COMPUTE n0 AGAIN
	#fourier_transform(lmda, lmda0, t)
	if(t_units == 'orbits'):
		#t *= 4.8658333 # to units of orbits (assuming constant orbital motion at 1751.7 degrees per day)
		t *= (n0 / 360) # to units of orbits (using the calculated mean orbital motion for this simulation) (days cancel, 360 deg = 1 orbit)
	else:
		t /= 365.25 # to units of years
	R_H = hill_radius(a0, masses[0])
	delta_a = compute_delta_a(a, a0, R_H, units)
	delta_lambda = compute_delta_lambda(lmda, lmda0)
	initial_sep_deg = compute_initial_separation_deg(lmda, lmda0)
	delta_e = compute_delta_e(e)
	#min_sep_deg1, min_sep_idx = compute_minimum_separation_deg(lmda, lmda0)
	#compute_deviation(lmda, lmda0)
	a_plot(bodies, a, a0, t, delta_a, R_H, units, t_units, black)
	E_plot(bodies, E, t, t_units, black)
	print('Mean energy: ' + str(np.mean(E)) + ' J')
	lmda_amp, min_sep_deg = summary_plot(bodies, t, a, a0, delta_a, R_H, e, delta_e, lmda, lmda0, delta_lambda, x, y, black, units, t_units)
	plot_longitude_histograms(lmda, False, black)
	record_data(masses, delta_a, lmda_amp, min_sep_deg, initial_sep_deg, np.mean(E), np.mean(L), units)
#main(black=False)
#longitude_histogram_from_file()