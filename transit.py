# Transit
#
# Author: Joseph A'Hearn
# Created 02/15/2021
#
# This program generates synthetic transit data based on a mercury6 simulation
#   

import numpy as np 
import bodies as bd 
import data_loader as dl 
import matplotlib.pyplot as plt 
import plot_assistant as pa 
import constants_of_mercury6 as cm6 
import random as rd 

def initialize_matrices(n_bodies, n_data):
	#return np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data))
	#return np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data))
	return np.zeros((n_bodies, n_data)), np.zeros((n_bodies, n_data))

def make_plot(t, f, t_min, t_max, i, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 't [days]', 'flux', t_min, t_max, np.amin(f) * 0.999, np.amax(f) * 1.001)
	ax.scatter(t, f, c='r')
	ax.plot([t_min, t_max], [1, 1], c='gray', linestyle='dashed')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='flux_vs_t_' + str(i) + '_' + str(int(np.round((t_max + t_min)/2))), black=black) 

def write_data(t, f):
	with open('transit_data.txt', 'a') as myfile:
	    myfile.write('t (days)\t') 
	    myfile.write('normalized flux\t\n')
	    for i in range(len(t)):
	    	myfile.write(str(t[i]) + ' \t\t')
	    	myfile.write(str(f[i]) + ' \n')

def main(filename='big.in', AU=cm6.AU()):
	bodies = bd.full_body_list(filename)
	masses = bd.full_mass_list(filename)
	# get the data from the aei files 
	t0 = dl.time_data(bodies[0])
	x, y = initialize_matrices(len(bodies), len(t0))

	data_surplus = 1
	#n_timesteps = 10
	
	t = np.linspace(t0[0], t0[-1], len(t0) * data_surplus)
	f = np.ones(len(t0) * data_surplus) # f is for flux; 1 is the default value, assuming that a transit is not taking place 
	# noise
	noise = 0.001
	for i in range(len(f)):
		f[i] += rd.uniform(-noise, noise)
	#x, y, u, v, lmda = initialize_matrices(len(bodies), len(t))
	#x, y, u, v, lmda, e, a, n = initialize_matrices(len(bodies), len(t))

	for i in range(len(bodies)):
		x_i, y_i = dl.xy_data(bodies[i])
		#x_i, y_i, u_i, v_i = dl.xyuv_data(bodies[i])
		for j in range(len(x_i)):
			x[i][j] = x_i[j] * AU
			y[i][j] = y_i[j] * AU
			#u[i][j] = u_i[j]
			#v[i][j] = v_i[j]

	M = masses[0]
	
	M_E = 5.9722E+24 # mass of Earth in kg
	
	if((0.1 * M_E) < M < (5 * M_E)):
		rho = 4000 # silicate
	elif((5 * M_E) < M < (25 * M_E)):
		rho = 1100 # ice giant
	elif(M > (25 * M_E)):
		rho = 900 # gas giant
	elif(M < (0.1 * M_E)):
		rho = 2600 # in between ice and silicate
	print(rho)
	R_planet = ((3 * M) / (4 * np.pi * rho))**(1/3)
	
	R_star = bd.R() # in meters 
	
	Delta_f = (R_planet / R_star)**2 # transit depth, equal to fraction of cross-sectional area covered
	print(Delta_f)
	
	Delta_y1 = R_star - R_planet # allowable distance from x-axis that still results in a total transit
	Delta_y2 = R_star + R_planet # allowable distance from x-axis that still results in a partial transit  

	print(Delta_y1, Delta_y2)
	
	for i in range(len(bodies)):
		for j in range(5, len(t0) - 6):
			if((x[i][j] > 0) and (np.sign(y[i][j]) != np.sign(y[i][j-1]))):
				#print('new transit')
				# go back to before the transit
				k = 0
				rewind = True
				while(rewind): 
					if(np.absolute(y[i][j+k]) < Delta_y2):
						if(np.absolute(y[i][j+k]) < Delta_y1): # complete transit
							f[j+k] -= Delta_f
						else: # partial transit
							#C = ((np.absolute(y[i][j+k]) - (R_star - R_planet)) / R_planet)**2 # should always be between 0 and 1; should correspond to the fraction of the area of the planet that is not covering the star
							C = ((np.absolute(y[i][j+k]) - Delta_y1) / (Delta_y2 - Delta_y1))**2 # this is an approximation
							f[j+k] -= Delta_f * (1 - C)
							
						k -= 1
					else:
						begidx = j + k - 5
						rewind = False
				
				# go forward to after the transit
				#print('going forward')
				k = 1
				ffw = True
				while(ffw):
					if(np.absolute(y[i][j+k]) < Delta_y2):
						if(np.absolute(y[i][j+k]) < Delta_y1): # complete transit
							f[j+k] -= Delta_f
						else: # partial transit
							#C = ((np.absolute(y[i][j+k]) - (R_star - R_planet)) / R_planet)**2 # should always be between 0 and 1; should correspond to the fraction of the area of the planet that is not covering the star
							C = ((np.absolute(y[i][j+k]) - Delta_y1) / (Delta_y2 - Delta_y1))**2 # this is an approximation
							f[j+k] -= Delta_f * (1 - C) 
							
						k += 1
					else:
						endidx = j + k + 6
						ffw = False
					

				#yk = np.linspace(y[i][j-1], y[i][j], data_surplus * n_timesteps)
				#begidx = int((j * data_surplus) - (data_surplus * n_timesteps / 2))
				#endidx = int((j * data_surplus) + (data_surplus * n_timesteps / 2))
				#for k in range(data_surplus * n_timesteps):
				#	if(i == 0):
				#		print(np.absolute(yk[k]) / Delta_y1)
				#	if(np.absolute(yk[k]) < Delta_y1): 
				#		#print('complete transit')
				#		f[(j * data_surplus) - data_surplus + k] -= Delta_f
				#		f[begidx + k] -= Delta_f
				#	elif(np.absolute(yk[k]) < Delta_y2): 
				#		#print('partial transit')
				#		C = ((np.absolute(yk[k]) - R_star) / R_planet)**2 # should always be between 0 and 1; should correspond to the fraction of the area of the planet that is not covering the star
				#		f[(j * data_surplus) - data_surplus + k] -= Delta_f * (1 - C) 
				#		f[begidx + k] -= Delta_f * (1 - C) 

				# later on, add occultations and noise, maybe even asteroseismology 

				
				make_plot(t[begidx:endidx], f[begidx:endidx], t[begidx], t[endidx], i)
	write_data(t, f)

#main()
#import sys 
#print("You are using Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))