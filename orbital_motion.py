# Orbital Motion
#
# Author: Joseph A'Hearn
# Created 06/26/2017
#
# This program includes functions that I could use often
#    when analyzing orbital motion
#   

import numpy as np 
import matplotlib.pyplot as plt 
from bodies import mu, full_body_list, R, J2
import data_loader as dl 
import plot_assistant as pa 
from useful import decapitalize
from constants_of_mercury6 import AU
#from xyz2geo import xyz2geo

def n_geo(body): # mean motion determined from the geometric orbital elements
	return 360 / (2 * np.pi * np.sqrt((dl.geo1mean(body))**3 / mu_Sat())) # deg/s

def n_rad(body): # mean motion determined from the geometric orbital elements
	return 1 / (np.sqrt((dl.geo1mean(body))**3 / mu_Sat())) # rad/s

def n_exact(mu1, mu2, body): # mean motion for each data point
	return 1 / (np.sqrt((dl.geo1data(body))**3 / (mu1 + mu2))) # rad/s

def T_geo(body): # mean motion determined from the geometric orbital elements
	return 2 * np.pi * np.sqrt((dl.geo1mean(body))**3 / mu_Sat())

def T_exact(mu1, mu2, body):
	return 2 * np.pi * np.sqrt((dl.geo1data(body))**3 / (mu1 + mu2))

def plot12(body):
	t = dl.geotime_data(body) / 365.25
	a, e = dl.geo12data(body)
	a *= 1.0E-03 # conversion to km
	fig, ax1, ax2 = pa.initialize21plot()
	fig, ax1 = pa.title_and_axes(fig, ax1, 'Semimajor axis', 't (yr)', r'$a$ (km)', t[0], t[-1], 0.999995 * np.amin(a), 1.000005 * np.amax(a))
	fig, ax2 = pa.title_and_axes(fig, ax2, 'Eccentricity',   't (yr)', r'$e$',      t[0], t[-1], 0, 1.08 * np.amax(e))
	ax1.scatter(t, a, s=0.1)
	ax2.scatter(t, e, s=0.1)
	plt.tight_layout()
	body = decapitalize(body)
	pa.save_and_clear_plot(fig, ax1, ax2, 'elements12_' + str(body))

def plot23(body):
	t = dl.time_data(body) / 365.25
	e, i = dl.geo23data(body)
	fig, ax1, ax2 = pa.initialize21plot()
	fig, ax1 = pa.title_and_axes(fig, ax1, 'Eccentricity', 't (yr)', 'e', t[0], t[-1], 0, 1.08 * np.amax(e))
	fig, ax2 = pa.title_and_axes(fig, ax2, 'Inclination',  't (yr)', 'i', t[0], t[-1], 0, 1.08 * np.amax(i))
	ax1.scatter(t, e, s=0.1)
	ax2.scatter(t, i, s=0.1)
	plt.tight_layout()
	body = decapitalize(body)
	pa.save_and_clear_plot(fig, ax1, ax2, 'elements23_' + str(body))

def plot56(body):
	t = dl.time_data(body) / 365.25
	o, n = dl.geo56data(body)
	fig, ax1, ax2 = pa.initialize21plot()
	fig, ax1 = pa.title_and_axes(fig, ax1, 'Longitude of pericenter',     't (yr)', r'$\varpi$', t[0], t[-1], 0, 360)
	fig, ax2 = pa.title_and_axes(fig, ax2, 'Longitude of ascending node', 't (yr)', r'$\Omega$', t[0], t[-1], 0, 360)
	ax1.scatter(t, o, s=0.1)
	ax2.scatter(t, n, s=0.1)
	ax1.set_yticks([0, 60, 120, 180, 240, 300, 360])
	ax2.set_yticks([0, 60, 120, 180, 240, 300, 360])
	plt.tight_layout()
	body = decapitalize(body)
	pa.save_and_clear_plot(fig, ax1, ax2, 'longitudes56_' + str(body))

def plot25(body):
	t = dl.time_data(body) / 365.25
	e, p = dl.geo25data(body)
	fig, ax1, ax2 = pa.initialize21plot()
	fig, ax1 = pa.title_and_axes(fig, ax1, 'Eccentricity', 't (yr)', 'e', t[0], t[-1], 0, 1.08 * np.amax(e))
	fig, ax2 = pa.title_and_axes(fig, ax2, 'Longitude of pericenter',  't (yr)', r'$\varpi$', t[0], t[-1], 0, 360)
	ax1.scatter(t, e, s=0.1)
	ax2.scatter(t, p, s=0.1)
	ax2.set_yticks([0, 60, 120, 180, 240, 300, 360])
	plt.tight_layout()
	body = decapitalize(body)
	pa.save_and_clear_plot(fig, ax1, ax2, 'elements25_' + str(body))

def plotgeo(body):
	t, a, e, i, m, p, l = dl.geo_data(body)
	t = t / 365.25
	fig, ax1, ax2, ax3, ax4, ax5, ax6 = pa.initialize32plot()
	fig, ax1 = pa.title_and_axes(fig, ax1, 'Semimajor axis', 't (yr)', 'a',                     t[0], t[-1], 0.995 * np.amin(a), 1.005 * np.amax(a))
	fig, ax2 = pa.title_and_axes(fig, ax2, 'Eccentricity',   't (yr)', 'e',                     t[0], t[-1], 0, 1.08 * np.amax(e))
	fig, ax3 = pa.title_and_axes(fig, ax3, 'Inclination',    't (yr)', 'i (deg)',               t[0], t[-1], 0, 1.08 * np.amax(i))
	fig, ax4 = pa.title_and_axes(fig, ax4, 'Mean longitude', 't (yr)', r'$\lambda$' + ' (deg)', t[0], t[-1], 0, 360)
	fig, ax5 = pa.title_and_axes(fig, ax5, 'Pericenter',     't (yr)', r'$\varpi$' + ' (deg)',  t[0], t[-1], 0.9 * np.amin(p), 1.1 * np.amax(p))
	fig, ax6 = pa.title_and_axes(fig, ax6, 'Ascending node', 't (yr)', r'$\Omega$' + ' (deg)',  t[0], t[-1], 0.9 * np.amin(l), 1.1 * np.amax(l))
	ax1.scatter(t, a, s=0.5)
	ax2.scatter(t, e, s=0.5)
	ax3.scatter(t, i, s=0.5)
	ax4.scatter(t, m, s=0.5)
	ax5.scatter(t, p, s=0.5)
	ax6.scatter(t, l, s=0.5)

	plt.tight_layout()
	body = decapitalize(body)
	pa.save_and_clear32plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, str(body) + '_geo_elements')

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

def define_corotating_frame(x, y):
	# theta is the angle at which the axes are rotated to keep body1 fixed
	theta = np.zeros(len(y))
	for j in range(len(y)):
		if(y[j]!= 0):
			theta[j] = (y[j] / np.absolute(y[j])) * np.arccos(x[j] / (np.sqrt(x[j]**2 + y[j]**2)))
		elif(x[j] < 0):
			theta[j] = np.pi 
		if(theta[j] < 0):
			theta[j] += 2 * np.pi 
	return theta

def convert_to_corotating_frame(theta, x, y): # this used to convert to Gm, but I removed that
	return ((x * np.cos(theta)) + (y * np.sin(theta))), ((-x * np.sin(theta)) + (y * np.cos(theta)))

def corotating_frame(fixed, free):
	x1, y1 = dl.xy_data(fixed) 
	x2, y2 = dl.xy_data(free) 
	theta = define_corotating_frame(x1, y1)
	x1c, y1c = convert_to_corotating_frame(theta, x1, y1)
	x2c, y2c = convert_to_corotating_frame(theta, x2, y2)
	#x2f = x2c - x1c # x position of body2 if we want body1 at the origin
	#x1f = 0 # if we want body1 at the origin
	return x1c, y1c, x2c, y2c

def plot_corotating_frame(fixed, free):
	x1c, y1c, x2c, y2c = corotating_frame(fixed, free)
	fig, ax = pa.initialize_plot()
	fig, ax = pa.title_and_axes(fig, ax, 'Corotating frame of ' + decapitalize(fixed), 'x (Gm)', 'y (Gm)')
	fig, ax = pa.Saturn_aerial(fig, ax)
	ax.scatter(x1c, y1c, s=0.1, label=decapitalize(fixed))
	ax.scatter(x2c, y2c, s=0.1, label=decapitalize(free))
	ax.legend(loc=2)
	pa.save_and_clear_plot(fig, ax, filename='corotating_frame_' + decapitalize(fixed)[0:3] + '_' + decapitalize(free)[0:3])

def plot_corotating_frame_N_equals_3(fixed, free1, free2, black=False, AU=AU()):
	x1c, y1c, x2c, y2c = corotating_frame(fixed, free1)
	x1c, y1c, x3c, y3c = corotating_frame(fixed, free2)
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, 'Corotating frame of ' + decapitalize(fixed), 'x [10' + r'$^3$' + ' km]', 'y [10' + r'$^3$' + ' km]')
	plt.gca().set_aspect('equal', adjustable='box')
	ax.scatter(x1c * AU / 1.0E+06, y1c * AU / 1.0E+06, s=0.1, label=decapitalize(fixed), color=assigned_color(2))
	ax.scatter(x2c * AU / 1.0E+06, y2c * AU / 1.0E+06, s=0.1, label=decapitalize(free1), color=assigned_color(3))
	ax.scatter(x3c * AU / 1.0E+06, y3c * AU / 1.0E+06, s=0.1, label=decapitalize(free2), color=assigned_color(1))
	ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='corotating_frame_' + decapitalize(fixed), black=black)

#plot_corotating_frame_N_equals_3('body0001','body0002', 'body0003')
#plot_corotating_frame_N_equals_3('body0002','body0003', 'body0001')
#plot_corotating_frame_N_equals_3('body0003','body0001', 'body0002')

def plot_positions_at_passes(fixed, passer, third_wheel):
	x1, y1 = dl.xy_data(fixed) 
	x2, y2 = dl.xy_data(passer)
	x3, y3 = dl.xy_data(third_wheel)
	theta = define_corotating_frame(x1, y1)
	x1c, y1c = convert_to_corotating_frame(theta, x1, y1)
	x2c, y2c = convert_to_corotating_frame(theta, x2, y2)
	x3c, y3c = convert_to_corotating_frame(theta, x3, y3)
	fig, ax = pa.initialize_plot()
	fig, ax = pa.title_and_axes(fig, ax, 'Corotating frame of ' + decapitalize(fixed), 'x (Gm)', 'y (Gm)', -0.3, 0.3, -0.3, 0.3)
	fig, ax = pa.Saturn_aerial(fig, ax)
	labeled = False
	for i in range(1, len(x1c)):
		if(x2c[i] > 0):
			if(np.sign(y2c[i]) != np.sign(y2c[i-1])): # then body1 and body2 just passed each other THIS IS ONLY RELIABLE IF THE STEP SIZE IS SMALL
				if(labeled == False):
					ax.scatter(x1c[i], y1c[i], s=0.1, c='r', label=decapitalize(fixed))
					ax.scatter(x2c[i], y2c[i], s=0.1, c='c', label=decapitalize(passer))
					ax.scatter(x3c[i], y3c[i], s=1, c='chartreuse', label=decapitalize(third_wheel))
					labeled = True
				else:	
					ax.scatter(x1c[i], y1c[i], s=0.1, c='r')
					ax.scatter(x2c[i], y2c[i], s=0.1, c='c')
					ax.scatter(x3c[i], y3c[i], s=0.1, c='chartreuse')
	ax.legend(loc=4)
	plt.gca().set_aspect('equal', adjustable='box')
	pa.save_and_clear_plot(fig, ax, filename='position_of_' + decapitalize(third_wheel)[0:3] + '_at_passes_of_' + decapitalize(fixed)[0:3] + '_and_'+ decapitalize(passer)[0:3])


#plot_positions_at_passes('PALLENE','MIMAS','ENCELADS')

def precession_rate_from_theory(body): # from equation 0.36 of Hedman's Introduction to Planetary Ring Dynamics
	return (3 / 2) * n_geo(body) * J2() * (R_Sat() / dl.geo1mean(body))**2

def plot_precession_rates(body):
	t = dl.geotime_data(body) / 365.25
	theory = precession_rate_from_theory(body) 
	print(body, theory)
	# the positive value is the time averaged rate of change of the longitude of pericenter
	# the negative value is the time averaged rate of change of the longitude of the ascending node
	o, n = dl.geo56data(body)
	odot = np.zeros(len(t))
	ndot = np.zeros(len(t))
	for i in range(1, len(o)):
		odot[i] = (o[i] - o[i-1]) / (t[i] - t[i-1]) # may have to insert modulo
		ndot[i] = (n[i] - n[i-1]) / (t[i] - t[i-1])
	print(np.mean(odot))
	print(np.mean(ndot))
	fig, ax1, ax2 = pa.initialize21plot()
	fig, ax1 = pa.title_and_axes(fig, ax1, 'Pericenter precession rate',     't (yr)', r'$\dot{\varpi}$', t[0], t[-1], -360, 360)
	fig, ax2 = pa.title_and_axes(fig, ax2, 'Ascending node precession rate', 't (yr)', r'$\dot{\Omega}$', t[0], t[-1], -360, 360)
	ax1.scatter(t, odot, s=0.5)
	ax2.scatter(t, ndot, s=0.5)
	#ax1.set_yticks([-60, -30, 0, 30, 60])
	#ax2.set_yticks([-60, -30, 0, 30, 60])
	plt.tight_layout()
	body = decapitalize(body)
	pa.save_and_clear_plot(fig, ax1, ax2, 'precession_rates56_' + str(body))


#body_list = full_body_list()
#for body in body_list:
#	plot_precession_rates(body)
#	plot12(body)
#	plot56(body) 



def true_anomaly(a, e, r, y): # y is a coordinate in the corotating frame with the pericenter as the fixed point
	# equation 2.139 in Solar System Dynamics
	cosf = (((a * (1 - (e**2))) / r) - 1) / e
	f = np.zeros(len(cosf))
	for i in range(len(cosf)):
		if(cosf[i] > 1):
				cosf[i] = 1
		if(cosf[i] < -1):
				cosf[i] = -1
		if(y[i] > 0):
			f[i] = np.arccos(cosf[i])
		else: 
			f[i] = (2 * np.pi) - np.arccos(cosf[i])
	return f # returns the value in radians

#def true_anomaly2(i, p, l, r, z): does not work
#	# i = inclination
#	# p = longitude of pericenter
#	# l = longitude of ascending node 
#	sinwf = z / (r * np.sin(np.deg2rad(i)))
#	f = np.zeros(len(sinwf))
#	for j in range(len(sinwf)):
#		if(sinwf[j] > 1):
#			sinwf[j] = 1
#		if(sinwf[j] < -1):
#			sinwf[j] = -1
#	return np.rad2deg(np.arcsin(sinwf)) + l - p # returns the value in degrees 

def plot_true_anomaly(body):
	a, e, i, p, n = dl.geo12356data(body)
	t, x, y, z, x_dot, y_dot, z_dot = dl.aei_data(body) 
	r = np.sqrt((x**2) + (y**2) + (z**2)) * AU()
	# convert to corotating frame with pericenter as the fixed point
	xc, yc = convert_to_corotating_frame(np.deg2rad(p), x, y)
	f = true_anomaly(a, e, r, yc)
	#f2 = true_anomaly2(i, p, n, r, z)
	fig, ax = pa.initialize_plot()
	fig, ax = pa.title_and_axes(fig, ax, 'True anomaly', 't (day)', 'f (deg)', t[0], t[-1], 0, 360)
	ax.scatter(t, f * 180 / np.pi, s=0.5, label='algorithm1')
	#ax.scatter(t, f2,              s=0.5, label='algorithm2')
	plt.tight_layout()
	plt.legend(loc=0)
	body = decapitalize(body)
	pa.save_and_clear_plot(fig, ax, filename='true_anomaly_' + str(body))

def true_anomaly_time_derivative(n, a, e, r, f): # from equation 2.32 in Solar System Dynamics
	return n * (a / r) * (1 + (e * np.cos(f))) / (np.sqrt(1 - (e**2)))

def true_anomaly_time_derivative2(n, a, e, r): # from the text above equation 2.31 in Solar System Dynamics
	return n * ((a / r)**2) * np.sqrt(1 - (e**2)) 

#plot_true_anomaly('MIMAS')
#plotgeo('AEGAEON')
#plotgeo('MIMAS')

def hkplot(body):
	e, p = dl.geo25data(body)
	fig, ax = pa.initialize_plot()
	fig, ax = pa.title_and_axes(fig, ax, 'Eccentricity Evolution', 'h', 'k', -np.amax(e), np.amax(e), -np.amax(e), np.amax(e))
	ax.scatter(e * np.cos(np.deg2rad(p)), e * np.sin(np.deg2rad(p)), s=0.5)
	pa.save_and_clear_plot(fig, ax, filename='eccentricity_evolution_' + str(decapitalize(body)))

def func1(e, ff, fi):
	return (1 / (1 + (e * np.cos(ff)))) - (1 / (1 + (e * np.cos(fi))))

def func2(e, ff, fi):
	return 0.5 * (np.sqrt(1 + (2 * e * np.cos(ff)) + (e**2)) + np.sqrt(1 + (2 * e * np.cos(fi)) + (e**2)))


def apoperiratio(e):
	# ratio of time spent near apoapsis to time spent near periapsis
	df = 0.5 * np.pi / 180 # half a degree in radians
	ffperi =  df 
	fiperi = 0
	ffapo  = np.pi + df
	fiapo  = np.pi 
	delta_t_peri = func1(e, ffperi, fiperi) / func2(e, ffperi, fiperi)
	delta_t_apo  = func1(e, ffapo,  fiapo ) / func2(e, ffapo,  fiapo )
	return np.absolute(delta_t_apo / delta_t_peri)

def plot_apoperiratio(emin=0.0001, emax=0.03, N=100, black=False):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, 'Ratio of time near apoapsis \nto time near periapsis', r'$e$', r'$\frac{T_{apo}}{T_{peri}}$', 0, emax, 1)
	e = np.linspace(emin, emax, N)
	ratio = np.zeros(N)
	for i in range(N):
		ratio[i] = apoperiratio(e[i])
	ax.plot(e, ratio)
	#ax.set_yscale('log')
	plt.tight_layout()
	pa.save_and_clear_plot(fig, ax, filename='apoperiratio_vs_e', black=black)




#print(apoperiratio(0.0196))
#plot_apoperiratio()
#plot12('AEGAEON')
#plot25('AEGAEON')
#plot25('4thbody')
