# D68 Resonance Search
#
# Author: Joseph A'Hearn
# Created 12/03/2018 for Pallene
# Modified 10/11/2019 for D68
#
# This program searches for a D68 resonance 

import numpy as np 
import data_loader as dl 
import bodies as bd 
#import xyz2geo as x2g 
import plot_assistant as pa 
import matplotlib.pyplot as plt 

def import_data(bodies, data_points):
	print("Importing geometric elements...")
	a = np.zeros((len(bodies), data_points))
	e = np.zeros((len(bodies), data_points))
	i = np.zeros((len(bodies), data_points))
	lmda = np.zeros((len(bodies), data_points))
	peri = np.zeros((len(bodies), data_points))
	node = np.zeros((len(bodies), data_points))
	for j in range(len(bodies)):
		aj, ej, ij, lmdaj, perij, nodej = dl.geo123456data(bodies[j])
		a[j] = aj
		e[j] = ej
		i[j] = ij
		lmda[j] = lmdaj
		peri[j] = perij
		node[j] = nodej
	return a, e, i, lmda, peri, node 

#def compute_frequencies(bodies, data_points, a, e, i, lmda, peri, node):
#	n = np.zeros((len(bodies), data_points))
#	kappa = np.zeros((len(bodies), data_points))
#	nu = np.zeros((len(bodies), data_points))
#	R_Sat = bd.R_Sat()
#	mu_Sat = bd.mu_Sat()
#	J2 = bd.J2()
#	J4 = bd.J4()
#	J6 = bd.J6()
#	for j in range(len(bodies)):
#		n[j] = x2g.mean_motion(R_Sat, mu_Sat, J2, J4, J6, (a[j], e[j], i[j], lmda[j], peri[j], node[j]))
#		kappa[j] = x2g.horizontal_epicyclic(R_Sat, mu_Sat, J2, J4, J6, (a[j], e[j], i[j], lmda[j], peri[j], node[j]))
#		nu[j] = x2g.vertical_epicyclic(R_Sat, mu_Sat, J2, J4, J6, (a[j], e[j], i[j], lmda[j], peri[j], node[j]))
#	return n, kappa, nu 

#def plot_resonant_argument(bodies, j, k, p, q, t, lmda, peri, node, res_type):
#	fig, ax = pa.initialize_plot()
#	fig, ax = pa.title_and_axes(fig, ax, '', r'$t$' + ' [yr]', r'$\phi$' + ' [deg]', t[0], t[-1], 0, 360)
#	if(res_type == 'e'):
#		phi = (((p + q) * lmda[j]) - (p * lmda[k]) - (q * peri[j])) % 360 
#		ax.set_title(str(p+q) + r'$\lambda_{' + (bodies[j]) + '}$' + 
#			' - ' +  str(p)   + r'$\lambda_{' + (bodies[k]) + '}$' + 
#			' - ' +  str(q)   + r'$\varpi_{'  + (bodies[j]) + '}$')
#	elif(res_type == 'i'):
#		phi = (((p + q) * lmda[j]) - (p * lmda[k]) - (q * node[j])) % 360
#		ax.set_title(str(p+q) + r'$\lambda_{' + (bodies[j]) + '}$' + 
#			' - ' +  str(p)   + r'$\lambda_{' + (bodies[k]) + '}$' + 
#			' - ' +  str(q)   + r'$\Omega_{'  + (bodies[j]) + '}$')
#	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
#	ax.scatter(t, phi, s=0.5)
#	plt.tight_layout()
#	pa.save_and_clear_plot(fig, ax, filename='resonant_argument_' + res_type + '_' + str(p+q) + bodies[j][0:3] + '_' + str(p) + bodies[k][0:3] + '_' + str(q) + bodies[j][0:3])

#def plot_resonant_argument2(bodies, j, k, m, p, q, s, u, w, t, lmda, peri, node):
#	fig, ax = pa.initialize_plot()
#	fig, ax = pa.title_and_axes(fig, ax, '', r'$t$' + ' [yr]', r'$\phi$' + ' [deg]', t[0], t[-1], 0, 360)
#	
#	phi = (((m) * lmda[j]) + (p * lmda[k]) + (q * peri[j]) + (s * peri[k]) + (u * node[j]) + (w * node[k])) % 360
#	ax.set_title(str(m) + r'$\lambda_{' + (bodies[j][0:1]) + '}$' + 
#		' + ' +  str(p) + r'$\lambda_{' + (bodies[k][0:1]) + '}$' + 
#		' + ' +  str(q) + r'$\varpi_{'  + (bodies[j][0:1]) + '}$\n' +
#		' + ' +  str(s) + r'$\varpi_{'  + (bodies[k][0:1]) + '}$' + 
#		' + ' +  str(u) + r'$\Omega_{'  + (bodies[j][0:1]) + '}$' + 
#		' + ' +  str(w) + r'$\Omega_{'  + (bodies[k][0:1]) + '}$' )
#	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
#	ax.scatter(t, phi, s=0.5)
#	plt.tight_layout()
#	pa.save_and_clear_plot(fig, ax, filename='resonant_argument_' + bodies[j][0:3] + '_' + bodies[k][0:3] + '_' + str(m) + '_' + str(p) + '_' + str(q) + '_' + str(s) + '_' + str(u) + '_' + str(w))

def plot_resonant_argument3(bodies0, bodies, j, k, l, m, p, b, q, s, d):
	# import data
	t = dl.geotime_data(bodies0[0])
	lmda_D68 = dl.geo4data(bodies0[0])
	peri_D68 = dl.geo5data(bodies0[0])
	if(l == 2): # Janus corresponds to l = 2	
		lmda_pert = dl.geo4data(bodies0[3])
		peri_pert = dl.geo5data(bodies0[3])
	elif(l == 3): # Mimas corresponds to l = 3
		lmda_pert = dl.geo4data(bodies0[1])
		peri_pert = dl.geo5data(bodies0[1])
		#m *= 7
		#p *= 7
		#b *= 7
		#s *= 7
		#d *= 7
	elif(l == 8): # Titan corresponds to l = 8
		lmda_pert = dl.geo4data(bodies0[8])
		peri_pert = dl.geo5data(bodies0[8])
		
	# Saturn
	lmda_Sat = 818.14 * t # time t is imported in days, so the rotation rate in degrees per day is easiest

	# Saturn offset?
	# D68 offset?
	t /= 365.25 # converting time from days to yrs




	phi = (((m) * lmda_Sat) + (p * lmda_D68) + (b * lmda_pert) + (s * peri_D68) + (d * peri_pert)) % 360

	fig, ax = pa.initialize_plot()
	fig, ax = pa.title_and_axes(fig, ax, '', r'$t$' + ' [yr]', r'$\phi$' + ' [deg]', t[0], t[-1], 0, 360)
	
	
	ax.set_title(str(m) + r'$\Omega_{' + (bodies[j][0:3]) + '}$' + 
		' + ' +  str(p) + r'$\lambda_{' + (bodies[k][0:3]) + '}$' + 
		' + ' +  str(b) + r'$\lambda_{' + (bodies[l][0:3]) + '}$\n' + 
		#' + ' +  str(q) + r'$\varpi_{'  + (bodies[j][0:1]) + '}$' +
		' + ' +  str(s) + r'$\varpi_{'  + (bodies[k][0:3]) + '}$' + 
		' + ' +  str(d) + r'$\varpi_{'  + (bodies[l][0:3]) + '}$')
	ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
	ax.scatter(t, phi, s=0.5)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='resonant_argument_' + bodies[j][0:3] + '_' + bodies[k][0:3] + '_' + bodies[l][0:3] + '_' + str(m) + '_' + str(p) + '_' + str(b) + '_' + str(q) + '_' + str(s) + '_' + str(d))

def general_two_body_resonance(n1, n2, pp1, pp2, np1, np2, m, p, q, s, u, w):
	return (m * n1) + (p * n2) + (q * pp1) + (s * pp2) + (u * np1) + (w * np2) 

def general_three_body_resonance(n1, n2, n3, pp1, pp2, pp3, np1, np2, np3, m, p, b, q, s, d, u, w, g):
	return (m * n1) + (p * n2) + (b * n3) + (q * pp1) + (s * pp2) + (d * pp3) + (u * np1) + (w * np2) + (g * np3)

def planar_three_body_resonance(n1, n2, n3, pp1, pp2, pp3, m, p, b, q, s, d):
	return (m * n1) + (p * n2) + (b * n3) + (q * pp1) + (s * pp2) + (d * pp3)

def mean_apsidal_precession(a, n, J2=0.01629071, R=60268):
	return 1.5 * n * J2 * ((R / a)**2) # averaged value over an orbit, from Hedman's Intro to Rings Chapter

def main(tol=1.0, maxrange=15, maxorder=4): # maxorder needs to be an even number for the coefficients of the nodal rates to be even numbers
	bodies0 = bd.full_body_list()
	bodies = ['Saturn', 'D68', 'Janus', 'Mimas', 'Enceladus', 'Tethys', 'Dione', 'Rhea', 'Titan']
	n = [818.14, 1735, 518.3431513, 381.9944948, 262.7318978, 190.6979109, 131.5349307, 79.6900459, 22.5769756] # units of degrees / day
	a = [60268, 67630, 151450, 185539, 238042, 294672, 377415, 527068, 1221865]
	omega_dot = np.zeros(len(n))
	kappa = np.zeros(len(n))
	for i in range(len(n)):
		omega_dot[i] = mean_apsidal_precession(a[i], n[i])
		kappa[i] = n[i] - omega_dot[i]
	#t = dl.geotime_data('PALLENE')
	#print('length of t: ' + str(len(t)))
	#a, e, i, lmda, peri, node = import_data(bodies, len(t))
	#n, kappa, nu = compute_frequencies(bodies, len(t), a, e, i, lmda, peri, node)
	#for j in range(len(bodies)):
	#	for k in range(len(bodies)):
	#		if(j != k):
	#			if((j == 1) or (k == 1)): 

	
	j = 0 # SATURN corresponds to j,k = 0
	k = 1 # D68 corresponds to j,k = 1
	for l in range(2, len(n)):
	# Janus corresponds to l = 2
	# Mimas corresponds to l = 3
	# Enceladus corresponds to l = 4
	# Tethys corresponds to l = 5
	# Dione corresponds to l = 6
	# Rhea corresponds to l = 7
	# Titan corresponds to l = 8

		#print('Corotation resonance search...')
		#for p in range(1, 20):
		#	for q in range(1, 5):
		#		# Two-body mean motion resonance
		#		phi_dot = ((p + q) * np.mean(n[l])) - (p * np.mean(n[k])) - (q * np.mean(kappa[l]))
		#		if(np.absolute(phi_dot) < tol):
		#			print(str(p + q) + ' n_' + bodies[l] + ' - ' + str(p) + ' n_' + bodies[k]  + ' - ' + str(q) + ' kappa_' + bodies[l] + ' = ' + str(phi_dot))
		#			#	plot_resonant_argument(bodies, j, k, p, q, t, lmda, peri, node, 'e')


	# TWO-BODY RESONANCE SEARCH
	#for m in range(1, maxrange):
	#	for p in range(-maxrange, 0):
	#		if((p != 0) and (np.absolute(np.absolute(p) - np.absolute(m)) < maxorder + 1)):
	#			for q in range(-maxorder, maxorder):
	#				for s in range(-maxorder, maxorder):
	#					for u in range(-maxorder, maxorder, 2):
	#						w = -(m + p + q + s + u)
	#						if(np.absolute(w) < maxorder) and (w % 2 == 0):
	#							phi_dot = general_two_body_resonance(np.mean(n[j]), np.mean(n[k]), np.mean(kappa[j]), np.mean(kappa[k]), np.mean(nu[j]), np.mean(nu[k]), m, p, q, s, u, w)
	#							if(np.absolute(phi_dot) < tol):
	#								print(str(bodies[j][0:3]) + ' ' + str(bodies[k][0:3]) + ' ' + str(m) + ' ' + str(p) + ' ' + str(q) + ' ' + str(s) + ' ' + str(u) + ' ' + str(w) + ' -> ' + str(phi_dot))
	#								plot_resonant_argument2(bodies, j, k, m, p, q, s, u, w, t, lmda, peri, node)

		print('Three-body resonance search...')
		# THREE-BODY RESONANCE SEARCH
		for m in range(1, maxrange):
			for p in range(-maxrange, maxrange):
				if(p != 0):
					for b in range(-maxrange, maxrange):
						if(b != 0):
							q = 0
							#for q in range(-maxorder, maxorder):
							for s in range(-maxorder, maxorder):
								d = -(m + p + b + q + s) 			
								if(np.absolute(d) < maxorder):
									phi_dot = planar_three_body_resonance(n[j], n[k], n[l], kappa[j], kappa[k], kappa[l], m, p, b, q, s, d)
									if(np.absolute(phi_dot) < tol):
										print(str(bodies[j][0:3]) + ' ' + str(bodies[k][0:3]) + ' ' + str(bodies[l][0:3]) + ' ' + str(m) + ' ' + str(p) + ' ' + str(b) + ' ' + str(q) + ' ' + str(s) + ' ' + str(d) + ' ' + ' -> ' + str(phi_dot))
										plot_resonant_argument3(bodies0, bodies, j, k, l, m, p, b, q, s, d)


main()