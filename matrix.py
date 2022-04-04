# Matrix
# Author: Joseph A'Hearn
# Created 04/10/2019
#
# This program contains functions that are useful
#   for dealing with matrices

import numpy as np


def initialize_matrix(N):
	return np.zeros((N, N))

def print_matrix(M):
	for i in range(len(M)):
		print(M[i])
	print('')

def column(M, c):
	return [M[i][c] for i in range(len(M))]

def row(M, r):
	return M[r][:]
  
def height(M):
	return len(M)
  
def width(M):
	return len(M[0])

def one_of_many_solutions(M, m, n):
	x = np.zeros(n)
	if(n == 5):
		x[4] = 1
		x[3] = 1
		x[2] = -((M[2][3] * x[3]) + M[2][4]) / M[2][2]
		x[1] = -((M[1][2] * x[2]) + (M[1][3] * x[3]) + M[1][4]) / M[1][1]
		x[0] = -((M[0][1] * x[1]) + (M[0][2] * x[2]) + (M[0][3] * x[3]) + M[0][4]) / M[0][0]
	if(n == 4):
		print('Multiple solutions')
		x[3] = 1
		x[2] = -((M[2][3] * x[3])) / M[2][2]
		x[1] = -((M[1][2] * x[2]) + (M[1][3] * x[3])) / M[1][1]
		x[0] = -((M[0][1] * x[1]) + (M[0][2] * x[2]) + (M[0][3] * x[3])) / M[0][0]
	return x

def one_solution(n, M):
	x = np.zeros(n)
	if(n == 5):
		x[4] = 1 
		x[3] = -M[3][4] / M[3][3]
		x[2] = -(M[2][4] + (M[2][3] * x[3])) / M[2][2]
		x[1] = -((M[1][2] * x[2]) + (M[1][3] * x[3]) + M[1][4]) / M[1][1]
		x[0] = -((M[0][1] * x[1]) + (M[0][2] * x[2]) + (M[0][3] * x[3]) + M[0][4]) / M[0][0]
		#print('solution')
		#print(x)
	elif(n == 3):
		x[2] = 1
		x[1] = -M[1][2] / M[1][1]
		x[0] = -((M[0][1] * x[1]) + M[0][2]) / M[0][0]
	elif(n == 4):
		# option 1
		#x[3] = 1
		#x[2] = -M[2][3] / M[2][2]
		#x[1] = -(M[1][3] + (M[1][2] * x[2])) / M[1][1]
		#x[0] = -((M[0][1] * x[1]) + (M[0][2] * x[2]) + M[0][3]) / M[0][0]
		# option 2
		#x[2] = 1 
		#x[3] = -M[2][2] / M[2][3]
		#x[1] = -(M[1][2] + (M[1][3] * x[3])) / M[1][1]
		#x[0] = -((M[0][1] * x[1]) + (M[0][3] * x[3]) + M[0][2]) / M[0][0]
		# option 3
		#x[1] = 1
		#x[3] = -M[1][1] / (M[1][3] - (M[1][2] * M[2][3] / M[2][2]))
		#x[2] = -(M[2][3] * x[3]) / M[2][2]
		#x[0] = -((M[0][2] * x[2]) + (M[0][3] * x[3]) + M[0][1]) / M[0][0]
		# option 4
		x[0] = 1
		x[3] = M[0][0] * M[1][1] / ((M[1][3] * (M[0][1] - M[1][1])) + ((M[2][3] / M[2][2]) * ((M[1][1] * M[0][2]) - (M[0][1] * M[1][2]))))
		x[1] = (x[3] / M[1][1]) * ((M[1][2] * M[2][3] / M[2][2]) - (M[1][3])) 
		x[2] = -(M[2][3] * x[3]) / M[2][2]
	else:
		x[n-1] = 1
		for i in range(n-1, -1, -1):
			s = 0
			for j in range(i, n - 1):
				s += M[i][j] * x[j]
			x[i] = (M[i][n - 1] - s) / M[i][i]
	return x 

def set_pivot(M, n, piv):
	M_new = initialize_matrix(n)
	for i in range(n):
		M_new[i] = M[piv[i]-1]
		#for j in range(n):
		#	if(piv[i] == j + 1):
		#		M_new[i] = M[i]
	return M_new, piv

def gaussian_elimination_with_set_pivot(M, piv, raiseValueError=False):
	n = height(M)
	M, m = set_pivot(M, n, piv)
	for i in range(n):
		for j in range(i+1, n):
			if(M[i][i] == 0): # then the matrix has at least two rows that are zeros
				#print('Multiple solutions')
				x = one_of_many_solutions(M, m, n)
				return x, m
			else:
				M[j] = [M[j][k] - M[i][k] * M[j][i] / M[i][i] for k in range(n)]
	print_matrix(M)
	print(m)
	if(M[n-1][n-1] == 0):
		if(raiseValueError == True):
			raise ValueError('No unique solution')
		else:
			x = one_solution(n, M)
			#print('No unique solution')
			# then we have 4 equations and 5 unknowns, but we set 1 unknown to 1, so then we have 4 equations and 4 unknowns
			
	elif(M[n-1][n-1] < 1.0E-05):
		x = one_solution(n, M)

	else:
		# backward substitution
		x = [0] * n
		for i in range(n-1, -1, -1):
			s = sum(M[i][j] * x[j] for j in range(i, n))
			x[i] = (M[i][n-1] - s) / M[i][i]
		#print('solution')
		#print(x)

	return x, m

def gaussian_elimination_with_pivot(M, m, raiseValueError=False):
	"""
	Parameters
	----------
	M  : list of list of floats (matrix)
	m  : list of ints or floats

	Returns
	-------  
	list of floats
	    solution to the system of linear equation
	list of ints or floats
	
	Raises
	------
	ValueError
	    no unique solution
	"""
	# forward elimination
	n = height(M)
	for i in range(n):
		pivot(M, n, i, m)
		for j in range(i+1, n):
			if(M[i][i] == 0): # then the matrix has at least two rows that are zeros
				#print('Multiple solutions')
				x = one_of_many_solutions(M, m, n)
				return x, m
			else:
				M[j] = [M[j][k] - M[i][k] * M[j][i] / M[i][i] for k in range(n)]
	print_matrix(M)
	print(m)
	
	if(M[n-1][n-1] == 0):
		if(raiseValueError == True):
			raise ValueError('No unique solution')
		else:
			x = one_solution(n, M)
			#print('No unique solution')
			# then we have 4 equations and 5 unknowns, but we set 1 unknown to 1, so then we have 4 equations and 4 unknowns
			
	elif(M[n-1][n-1] < 1.0E-05):
		x = one_solution(n, M)

	else:
		# backward substitution
		x = [0] * n
		for i in range(n-1, -1, -1):
			s = sum(M[i][j] * x[j] for j in range(i, n))
			x[i] = (M[i][n-1] - s) / M[i][i]
		#print('solution')
		#print(x)

	return x, m
'''
# shorter way to pivot but cannot run in trinket
def pivot(m, n, i):
  max_row = max(range(i, n), key=lambda r: abs(m[r][i]))
  m[i], m[max_row] = m[max_row], m[i]
'''

def pivot(M, n, i, m):
	pivoted = False
	maximum = -1e100
	for r in range(i, n):
		if(maximum < abs(M[r][i])):
			max_row = r
			maximum = abs(M[r][i])
			pivoted = True
	if(pivoted == True):
		M[i], M[max_row] = M[max_row], M[i]
		m[i], m[max_row] = m[max_row], m[i]


#def row_reduction_1(M):
#	if(len(M[0]) == 5):
#		M1 = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
#		M1[0] = M[0]
#		M1[1] = [M[1][0] - ((M[1][0]/M[0][0]) * M[0][0]), M[1][1] - ((M[1][0]/M[0][0]) * M[0][1]), M[1][2] - ((M[1][0]/M[0][0]) * M[0][2]), M[1][3] - ((M[1][0]/M[0][0]) * M[0][3]), M[1][4] - ((M[1][0]/M[0][0]) * M[0][4])]
#		M1[2] = [M[2][0] - ((M[2][0]/M[0][0]) * M[0][0]), M[2][1] - ((M[2][0]/M[0][0]) * M[0][1]), M[2][2] - ((M[2][0]/M[0][0]) * M[0][2]), M[2][3] - ((M[2][0]/M[0][0]) * M[0][3]), M[2][4] - ((M[2][0]/M[0][0]) * M[0][4])]
#		M1[3] = [M[3][0] - ((M[3][0]/M[0][0]) * M[0][0]), M[3][1] - ((M[3][0]/M[0][0]) * M[0][1]), M[3][2] - ((M[3][0]/M[0][0]) * M[0][2]), M[3][3] - ((M[3][0]/M[0][0]) * M[0][3]), M[3][4] - ((M[3][0]/M[0][0]) * M[0][4])]
#		M1[4] = [M[4][0] - ((M[4][0]/M[0][0]) * M[0][0]), M[4][1] - ((M[4][0]/M[0][0]) * M[0][1]), M[4][2] - ((M[4][0]/M[0][0]) * M[0][2]), M[4][3] - ((M[4][0]/M[0][0]) * M[0][3]), M[4][4] - ((M[4][0]/M[0][0]) * M[0][4])]
#	
#	# if a row already has a zero there, then we don't need to do the row reduction for that row
#	return M1

def singularity_check(M):
	if(np.absolute(np.linalg.det(M)) < np.finfo(float).eps):
		print('The determinant of matrix M is ' + str(np.linalg.det(M)) + ', the absolute value of which is less than machine precision.')
		print('Thus the matrix is singular.')
	else:
		print('The determinant of matrix M is ' + str(np.linalg.det(M)))

def diagonal_check(M):
	need_to_pivot = False
	for i in range(len(M)):
		if (M[i][i] != 0):
			continue
		else:
			#print('This matrix must be pivoted.')
			need_to_pivot = True
			break
	return need_to_pivot

def reorder(x, m):
	y = np.zeros(len(m))
	for i in range(len(m)):
		for j in range(len(m)):
			if(m[i] == j + 1):
				y[j] = x[i]
	return y

