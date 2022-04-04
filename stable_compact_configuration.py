# Stable Compact Configuration
#
# Author: Joseph A'Hearn
# Created 05/24/2019
#
# This program reproduces Figure 3 from Salo and Yoder 1988
#   

import numpy as np 
import plot_assistant as pa 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib.path as path

def angles_of_configIa(N): # copied from ~Downloads/procyon/initial_positions.py
	if(N == 2):
		angles = [60.0]
	elif(N == 3):
		angles = [47.361]
	elif(N == 4):
		angles = [37.356, 41.498]
	elif(N == 5):
		angles = [32.660, 38.202]
	elif(N == 6):
		angles = [28.536, 30.013, 36.310]
	elif(N == 7):
		angles = [26.278, 28.535, 35.463]
	elif(N == 8):
		angles = [24.460, 25.248, 28.114, 35.902]
	else:
		print('Type Ia stationary configurations do not exist for N = ' + str(N))
		angles = 0
	return angles

def compute_halfspread(i, angles):
	if((i % 2) == 0):
		if(i == 2):
			halfspread = angles[0] / 2
		else:
			halfspread = np.sum(angles[1:]) + angles[0] / 2
	else:
		halfspread = np.sum(angles)
	return halfspread

def Nequals3(black=True, dpi=300, fs=16):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', '', '', -1, 1, -0.25, 1)
	plt.gca().set_aspect('equal', adjustable='box')
	if(black == True):
		c = 'w'
		o = 'k'
	else:
		c = 'k'
		o = 'w'
	ax.plot([0,0], [-0.05, 0.05], color=c)
	ax.plot([-0.05, 0.05], [0,0], color=c)
	i = 3
	angles = angles_of_configIa(i)
	r = 0.12 * i
	r_t = r + 0.0325
	halfspread = compute_halfspread(i, angles)
	arc = patches.Arc((0,0), 2 * r, 2 * r, angle=0.0, theta1=90-halfspread, theta2=90+halfspread, color=c, linestyle='-', linewidth=2)
	ax.add_artist(arc)
	angle = angles[0] / 2
	ax.scatter(0, r, color=c, s=50, zorder=3)
	ax.scatter(0, r, color=o, s=30, zorder=4)
	ax.scatter(0, r, color=c, s=10, zorder=5)
	for j in range(len(angles)):
		ax.scatter(r * np.cos(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=c, s=10, zorder=5)
		if(j != 0):
			angle += angles[j] 
		ax.text(r_t * np.cos(-np.deg2rad(angle) + (np.pi / 2)), r_t * np.sin(np.deg2rad(angle) + (np.pi / 2)), str(np.round(angles[j], 1)) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-angle, fontsize=fs)	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	if(black==True):
		#ax.spines['bottom'].set_color('white')
		#ax.spines['left'].set_color('white')
		#ax.tick_params(axis='x', colors='white')
		#ax.tick_params(axis='y', colors='white')
		#ax.tick_params(axis='both', direction='in')
		#ax.get_xaxis().tick_bottom()
		#ax.get_yaxis().tick_left()
		ax.figure.savefig("type_Ia_Nequals3_configuration_" + str(black) + ".png", dpi=dpi, facecolor=fig.get_facecolor(), transparent=True)
	else:
		ax.figure.savefig("type_Ia_Nequals3_configuration_" + str(black) + ".png", dpi=dpi)
	plt.clf()

def main(black=True, dpi=300, fs=16):
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', '', '', -1, 1, -0.25, 1)
	plt.gca().set_aspect('equal', adjustable='box')
	if(black == True):
		c = 'w'
		o = 'k'
	else:
		c = 'k'
		o = 'w'
	ax.plot([0,0], [-0.05, 0.05], color=c)
	ax.plot([-0.05, 0.05], [0,0], color=c)
	for i in range(2, 9):
		angles = angles_of_configIa(i)
		r = 0.12 * i
		r_t = r + 0.0325
		halfspread = compute_halfspread(i, angles)
		arc = patches.Arc((0,0), 2 * r, 2 * r, angle=0.0, theta1=90-halfspread, theta2=90+halfspread, color=c, linestyle='-', linewidth=2)
		ax.add_artist(arc)
		angle = angles[0] / 2
		if((i % 2) == 0): # 2,4,6,8
			for j in range(len(angles)):
				if(j != 0):
					angle += angles[j]
				ax.scatter(r * np.cos(+np.deg2rad(angle) + (np.pi / 2)), r * np.sin(+np.deg2rad(angle) + (np.pi / 2)), color=c, s=50, zorder=3)
				ax.scatter(r * np.cos(+np.deg2rad(angle) + (np.pi / 2)), r * np.sin(+np.deg2rad(angle) + (np.pi / 2)), color=o, s=30, zorder=4)
				ax.scatter(r * np.cos(+np.deg2rad(angle) + (np.pi / 2)), r * np.sin(+np.deg2rad(angle) + (np.pi / 2)), color=c, s=10, zorder=5)
				ax.scatter(r * np.cos(-np.deg2rad(angle) + (np.pi / 2)), r * np.sin(-np.deg2rad(angle) + (np.pi / 2)), color=c, s=50, zorder=3)
				ax.scatter(r * np.cos(-np.deg2rad(angle) + (np.pi / 2)), r * np.sin(-np.deg2rad(angle) + (np.pi / 2)), color=o, s=30, zorder=4)
				ax.scatter(r * np.cos(-np.deg2rad(angle) + (np.pi / 2)), r * np.sin(-np.deg2rad(angle) + (np.pi / 2)), color=c, s=10, zorder=5)
				if(j == 0):
					if(i == 2):
						ax.text(0, r_t, str(int(angles[j])) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, fontsize=fs)
					else:
						ax.text(0, r_t, str(np.round(angles[j], 1)) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, fontsize=fs)
				else:
					#if(j == 1):
					#	angulo = angles[0] + angles[1] / 2
					#elif(j == 2):
					#	angulo = angles[0] + angles[1] + (angles[2] / 2)
					#elif(j == 3):
					#	angulo = angles[0] + angles[1] + angles[2] + (angles[3] / 2)
					#ax.text(r_t * np.cos(-angulo + (np.pi / 2)), r_t * np.sin(-angulo + (np.pi / 2)), str(np.round(angles[j], 1)) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-np.sum(angles[:j]))
					ax.text(r_t * np.cos(-np.deg2rad(np.sum(angles[:j])) + (np.pi / 2)), r_t * np.sin(np.deg2rad(np.sum(angles[:j])) + (np.pi / 2)), str(np.round(angles[j], 1)) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-np.sum(angles[:j]), fontsize=fs)
			
		else: #3,5,7
			ax.scatter(0, r, color=c, s=50, zorder=3)
			ax.scatter(0, r, color=o, s=30, zorder=4)
			ax.scatter(0, r, color=c, s=10, zorder=5)
			for j in range(len(angles)):
				ax.scatter(r * np.cos(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=c, s=50, zorder=3)
				ax.scatter(r * np.cos(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=o, s=30, zorder=4)
				ax.scatter(r * np.cos(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(+np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=c, s=10, zorder=5)
				ax.scatter(r * np.cos(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=c, s=50, zorder=3)
				ax.scatter(r * np.cos(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=o, s=30, zorder=4)
				ax.scatter(r * np.cos(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), r * np.sin(-np.deg2rad(np.sum(angles[:j+1])) + (np.pi / 2)), color=c, s=10, zorder=5)
				if(j != 0):
					angle += angles[j] 
				ax.text(r_t * np.cos(-np.deg2rad(angle) + (np.pi / 2)), r_t * np.sin(np.deg2rad(angle) + (np.pi / 2)), str(np.round(angles[j], 1)) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-angle, fontsize=fs)	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	if(black==True):
		#ax.spines['bottom'].set_color('white')
		#ax.spines['left'].set_color('white')
		#ax.tick_params(axis='x', colors='white')
		#ax.tick_params(axis='y', colors='white')
		#ax.tick_params(axis='both', direction='in')
		#ax.get_xaxis().tick_bottom()
		#ax.get_yaxis().tick_left()
		ax.figure.savefig("type_Ia_configurations_" + str(black) + ".png", dpi=dpi, facecolor=fig.get_facecolor(), transparent=True)
	else:
		ax.figure.savefig("type_Ia_configurations_" + str(black) + ".png", dpi=dpi)
	plt.clf()

def D68(black=True, dpi=300, bw=False):
	fig, ax = pa.initialize_plot(black)
	#fig, ax = pa.title_and_axes(fig, ax, '', '', '', -0.5, 0.5, -0.1, 0.5)
	fig, ax = pa.title_and_axes(fig, ax, '', '', '', -1, 1, -0.25, 1)
	plt.gca().set_aspect('equal', adjustable='box')
	if(black == True):
		c = 'w'
		o = 'k'
	else:
		c = 'k'
		o = 'w'
	ax.plot([0,0], [-0.05, 0.05], color=c)
	ax.plot([-0.05, 0.05], [0,0], color=c)
	angles = [26, 32, 29]
	r = 0.12 * 4
	r_t = r + 0.0325
	offset = 11.5
	halfspread = np.sum(angles) / 2
	arc = patches.Arc((0,0), 2 * r, 2 * r, angle=0.0, theta1=90-halfspread, theta2=90+halfspread, color=c, linestyle='-', linewidth=2)
	ax.add_artist(arc)
	if(bw == True):
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
	else:
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color='m', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color='m', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color='c', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color='c', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color='g', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color='g', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color='darkorange', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color='darkorange', s=10, zorder=5)
	#ax.text(r_t *  np.cos(np.deg2rad(-43 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-43 + offset) + (np.pi / 2)), str(angles[0]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-30)
	#ax.text(r_t *  np.cos(np.deg2rad(-14 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-14 + offset) + (np.pi / 2)), str(angles[1]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=0)
	#ax.text(r_t *  np.cos(np.deg2rad( 16 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad( 16 + offset) + (np.pi / 2)), str(angles[2]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=30)
	ax.text(r_t *  np.cos(np.deg2rad(-43 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-43 + offset) + (np.pi / 2)), str(angles[0]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-30, fontsize=16)
	ax.text(r_t *  np.cos(np.deg2rad(-14 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-14 + offset) + (np.pi / 2)), str(angles[1]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=0  , fontsize=16)
	ax.text(r_t *  np.cos(np.deg2rad( 16 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad( 16 + offset) + (np.pi / 2)), str(angles[2]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=30 , fontsize=16)
	#ax.text(-0.3, 0.31, 'LL', horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text(-0.1, 0.43, 'L',  horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text( 0.12, 0.43, 'M',  horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text( 0.3, 0.31, 'T',  horizontalalignment='center', verticalalignment='center', color=c)
	# [-55, -29, 3, 32
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	if(black==True):
		#ax.spines['bottom'].set_color('white')
		#ax.spines['left'].set_color('white')
		#ax.tick_params(axis='x', colors='white')
		#ax.tick_params(axis='y', colors='white')
		#ax.tick_params(axis='both', direction='in')
		#ax.get_xaxis().tick_bottom()
		#ax.get_yaxis().tick_left()
		ax.figure.savefig("D68_configuration_" + str(black) + ".png", dpi=dpi, facecolor=fig.get_facecolor(), transparent=True)
	else:
		ax.figure.savefig("D68_configuration_" + str(black) + ".png", dpi=dpi)
	plt.clf()

def D68plus5(black=True, dpi=300, bw=False):
	fig, ax = pa.initialize_plot(black)
	#fig, ax = pa.title_and_axes(fig, ax, '', '', '', -0.5, 0.5, -0.1, 0.5)
	fig, ax = pa.title_and_axes(fig, ax, '', '', '', -1, 1, -0.25, 1)
	plt.gca().set_aspect('equal', adjustable='box')
	if(black == True):
		c = 'w'
		o = 'k'
		s = 'chartreuse'
	else:
		c = 'k'
		o = 'w'
		s = 'r'
	ax.plot([0,0], [-0.05, 0.05], color=c)
	ax.plot([-0.05, 0.05], [0,0], color=c)
	angles  = [26, 32, 29]
	angles2 = [32, 26, 32, 29, 33.5]
	r = 0.12 * 4
	r_t = r + 0.0325
	offset = 11.5
	halfspread  = np.sum(angles) / 2
	halfspread2 = np.sum(angles2) / 2
	arc  = patches.Arc((0,0), 2 * r, 2 * r, angle=0.0, theta1=90-halfspread,  theta2=90+halfspread,  color=c, linestyle='-', linewidth=2, zorder=2)
	arc2 = patches.Arc((0,0), 2 * r, 2 * r, angle=0.0, theta1=90-halfspread2, theta2=90+halfspread2, color=s, linestyle='--', linewidth=1, zorder=1)
	ax.add_artist(arc)
	ax.add_artist(arc2)
	if(bw == True):
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
	else:
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color='m', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(-55 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-55 + offset) + (np.pi / 2)), color='m', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color='c', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(-29 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(-29 + offset) + (np.pi / 2)), color='c', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color='g', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(  3 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(  3 + offset) + (np.pi / 2)), color='g', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color='darkorange', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad( 32 + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad( 32 + offset) + (np.pi / 2)), color='darkorange', s=10, zorder=5)
	#ax.text(r_t *  np.cos(np.deg2rad(-43 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-43 + offset) + (np.pi / 2)), str(angles[0]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-30)
	#ax.text(r_t *  np.cos(np.deg2rad(-14 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-14 + offset) + (np.pi / 2)), str(angles[1]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=0)
	#ax.text(r_t *  np.cos(np.deg2rad( 16 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad( 16 + offset) + (np.pi / 2)), str(angles[2]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=30)
	ax.text(r_t *  np.cos(np.deg2rad(-43 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-43 + offset) + (np.pi / 2)), str(angles[0]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-30, fontsize=16)
	ax.text(r_t *  np.cos(np.deg2rad(-14 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-14 + offset) + (np.pi / 2)), str(angles[1]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=0  , fontsize=16)
	ax.text(r_t *  np.cos(np.deg2rad( 16 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad( 16 + offset) + (np.pi / 2)), str(angles[2]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=30 , fontsize=16)
	#ax.text(-0.3, 0.31, 'LL', horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text(-0.1, 0.43, 'L',  horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text( 0.12, 0.43, 'M',  horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text( 0.3, 0.31, 'T',  horizontalalignment='center', verticalalignment='center', color=c)
	for i in np.arange(-96, -78, 0.5):
		ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color=s, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color=s, s=10, zorder=5)
	for i in np.arange(55, 76, 0.5):
		ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color=s, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color=s, s=10, zorder=5)
	ax.text( 0.5, 0.2, r'$32^{\circ} \pm 11^{\circ}$', horizontalalignment='center', verticalalignment='center', color=s, rotation=-65, fontsize=16)
	ax.text(-0.5, 0.2, r'$33^{\circ} \pm 11^{\circ}$', horizontalalignment='center', verticalalignment='center', color=s, rotation= 65, fontsize=16)

	# labels for the clumps
	ax.text(-0.31, 0.29, 'LL', color='darkorange', fontsize=16)
	ax.text(-0.12, 0.36, 'L', color='g', fontsize=16)
	ax.text(0.1, 0.36, 'M', color='c', fontsize=16)
	ax.text(0.25, 0.29, 'T', color='m', fontsize=16)


	# [-55, -29, 3, 32
	ax.axes.get_xaxis().set_visible(False)
	ax.axes.get_yaxis().set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	if(black==True):
		#ax.spines['bottom'].set_color('white')
		#ax.spines['left'].set_color('white')
		#ax.tick_params(axis='x', colors='white')
		#ax.tick_params(axis='y', colors='white')
		#ax.tick_params(axis='both', direction='in')
		#ax.get_xaxis().tick_bottom()
		#ax.get_yaxis().tick_left()
		ax.figure.savefig("D68_configuration_plus_5_" + str(black) + ".png", dpi=dpi, facecolor=fig.get_facecolor(), transparent=True)
	else:
		ax.figure.savefig("D68_configuration_plus_5_" + str(black) + ".png", dpi=dpi)
	plt.clf()

def D68plus5specific(black=True, dpi=300, bw=False):
	fig, ax = pa.initialize_plot(black)
	#fig, ax = pa.title_and_axes(fig, ax, '', '', '', -0.5, 0.5, -0.1, 0.5)
	fig, ax = pa.title_and_axes(fig, ax, '', '', '', -1, 1, -0.25, 1)
	plt.gca().set_aspect('equal', adjustable='box')
	if(black == True):
		c = 'w'
		o = 'k'
	else:
		c = 'k'
		o = 'w'
	ax.plot([0,0], [-0.05, 0.05], color=c)
	ax.plot([-0.05, 0.05], [0,0], color=c)
	#angles  = [26, 32, 29, 33] #no perturbation
	angles  = [21, 27, 24, 28]
	#angles2 = [32, 26, 32, 29, 33.5]
	r = 0.12 * 4
	r_t = r + 0.0325
	offset = -5

	# no perturbation
	#lmda_t  = -55
	#lmda_m  = -29
	#lmda_l  =   3
	#lmda_ll =  32
	#lmda_5  =  65

	lmda_t  = -45
	lmda_m  = -24
	lmda_l  =   3
	lmda_ll =  27
	lmda_5  =  55



	halfspread  = np.sum(angles) / 2
	#halfspread2 = np.sum(angles2) / 2
	arc  = patches.Arc((0,0), 2 * r, 2 * r, angle=0.0, theta1=90-halfspread,  theta2=90+halfspread,  color=c, linestyle='-', linewidth=2, zorder=2)
	#arc2 = patches.Arc((0,0), 2 * r, 2 * r, angle=0.0, theta1=90-halfspread2, theta2=90+halfspread2, color='chartreuse', linestyle='--', linewidth=1, zorder=1)
	ax.add_artist(arc)
	#ax.add_artist(arc2)
	if(bw == True):
		ax.scatter(r * np.cos(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), color=c, s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), color=c, s=10, zorder=5)
	else:
		ax.scatter(r * np.cos(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), color='m', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_t  + offset) + (np.pi / 2)), color='m', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), color='c', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_m  + offset) + (np.pi / 2)), color='c', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), color='g', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_l  + offset) + (np.pi / 2)), color='g', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), color='darkorange', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_ll + offset) + (np.pi / 2)), color='darkorange', s=10, zorder=5)
		ax.scatter(r * np.cos(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), color='r', s=50, zorder=3)
		ax.scatter(r * np.cos(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
		ax.scatter(r * np.cos(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(lmda_5  + offset) + (np.pi / 2)), color='r', s=10, zorder=5)
	#ax.text(r_t *  np.cos(np.deg2rad(-43 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-43 + offset) + (np.pi / 2)), str(angles[0]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=-30)
	#ax.text(r_t *  np.cos(np.deg2rad(-14 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(-14 + offset) + (np.pi / 2)), str(angles[1]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=0)
	#ax.text(r_t *  np.cos(np.deg2rad( 16 + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad( 16 + offset) + (np.pi / 2)), str(angles[2]) + r'$^{\circ}$' + ' ' + r'$\pm$' + ' ' r'$1^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=30)
	ax.text(r_t *  np.cos(np.deg2rad(((lmda_t + lmda_m ) / 2) + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(((lmda_t + lmda_m ) / 2) + offset) + (np.pi / 2)), str(angles[0]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=((lmda_t + lmda_m ) / 2), fontsize=16)
	ax.text(r_t *  np.cos(np.deg2rad(((lmda_m + lmda_l ) / 2) + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(((lmda_m + lmda_l ) / 2) + offset) + (np.pi / 2)), str(angles[1]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=((lmda_l + lmda_m ) / 2), fontsize=16)
	ax.text(r_t *  np.cos(np.deg2rad(((lmda_l + lmda_ll) / 2) + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(((lmda_l + lmda_ll) / 2) + offset) + (np.pi / 2)), str(angles[2]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=((lmda_l + lmda_ll) / 2), fontsize=16)
	ax.text(r_t *  np.cos(np.deg2rad(((lmda_5 + lmda_ll) / 2) + offset) + (np.pi / 2)), r_t * np.sin(np.deg2rad(((lmda_5 + lmda_ll) / 2) + offset) + (np.pi / 2)), str(angles[3]) + r'$^{\circ}$', horizontalalignment='center', verticalalignment='center', color=c, rotation=((lmda_5 + lmda_ll) / 2), fontsize=16)
	#ax.text(-0.3, 0.31, 'LL', horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text(-0.1, 0.43, 'L',  horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text( 0.12, 0.43, 'M',  horizontalalignment='center', verticalalignment='center', color=c)
	#ax.text( 0.3, 0.31, 'T',  horizontalalignment='center', verticalalignment='center', color=c)
	#for i in np.arange(-96, -78, 0.5):
	#	ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color='chartreuse', s=50, zorder=3)
	#	ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
	#	ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color='chartreuse', s=10, zorder=5)
	#for i in np.arange(55, 76, 0.5):
	#	ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color='chartreuse', s=50, zorder=3)
	#	ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color=o, s=30, zorder=4)
	#	ax.scatter(r * np.cos(np.deg2rad(i + offset) + (np.pi / 2)), r *   np.sin(np.deg2rad(i + offset) + (np.pi / 2)), color='chartreuse', s=10, zorder=5)
	#ax.text( 0.5, 0.2, r'$32^{\circ} \pm 11^{\circ}$', horizontalalignment='center', verticalalignment='center', color='chartreuse', rotation=-65)
	#ax.text(-0.5, 0.2, r'$33^{\circ} \pm 11^{\circ}$', horizontalalignment='center', verticalalignment='center', color='chartreuse', rotation= 65)

	# [-55, -29, 3, 32
	if(black==True):
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
		#ax.spines['bottom'].set_color('white')
		#ax.spines['left'].set_color('white')
		#ax.tick_params(axis='x', colors='white')
		#ax.tick_params(axis='y', colors='white')
		#ax.tick_params(axis='both', direction='in')
		#ax.get_xaxis().tick_bottom()
		#ax.get_yaxis().tick_left()
		ax.figure.savefig("D68_configuration_plus_5_specific_" + str(black) + ".png", dpi=dpi, facecolor=fig.get_facecolor(), transparent=True)
	else:
		ax.figure.savefig("D68_configuration_plus_5_specific_" + str(black) + ".png", dpi=dpi)
	plt.clf()

#main(black=False)
#D68(black=False)
#D68plus5(black=False)
#D68plus5specific(black=False)
Nequals3(False)
