import time
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import sympy as sym 
from scipy.optimize import curve_fit

def data_input(file_location):
	# Inputs data frame
	df = pd.read_csv(file)
	return df

def f (x, p0):
	# The equation to be fitted, however with the parameters stored in p0
	c, b1, m1, a, b2, m2 = p0
	return (c / (1 + np.exp(-b1 * (x - m1)))) * (a + (c * np.exp((-1 * (np.exp(-b2 * (x - m2)))))))

def f_fit (x, c, b1, m1, a, b2, m2):
	#The equation to be fitted using the curve_fit protocol, in contrast with f, the parameters
	# are not packed, and are available in a form that curve_fit can manipulate
	return (c / (1 + np.exp(-b1 * (x - m1)))) * (a + (c * np.exp((-1 * (np.exp(-b2 * (x - m2)))))))

def residual_calculator(tt, y, p0):
	#Error checking, make sure input is not a string
	if type(p0[0]) == str:
		y_est = [x for x in range(0, 21)]
		res2 = [x for x  in range(0, 21)]
		return y_est, res2
	
	#This part of the function will calculate the residual
	y_est = []
	res2 = []
	for a in tt:
		count = 0
		#Calculates the function with the given parameters
		temp = f(a, p0)
		y_est.append(temp)
		temp_res = np.power((y[count] - temp), 2)
		res2.append(temp_res)
		count += 1
	return y_est, res2

def marquardt(f, x, y, p0):
	# This module fits the data points to the curve using 
	# optimize.curve_fit from scipy
	try:
		p, cov = curve_fit(f, x, y, p0, maxfev=1000)
		out = [p[0], p[1], p[2], p[3], p[4], p[5]]
	except RuntimeError:
		#If curve_fit can't fit the equation to the data points
		out = ["NaN", "NaN", "NaN", "NaN", "NaN", "NaN"]
	return out

def plot_rfr(p0):
	if type(p0[0]) == str:
		y_pred = 'NaN'
		return y_pred
	y_pred = [f(a, p0) for a in xrange(0, 3001)]
	return y_pred

def max_rfr(rfr_list):
	# Finds the maximum r:fr and the tt to the max r:fr value 
	max_rfr = max(rfr_list)
	max_rfr_tt = rfr_list.index(max_rfr)
	return max_rfr, max_rfr_tt

def senescence(max_rfr_tt, sen_rfr, rfr_list):
	'''Perform type check before this'''
	# Calculates the tt to senescence
	for item in (rfr_list[(rfr_list.index(max)):]):
		if round(item, 2) == round(item, 2):
			sen = item
			break
		else:
			sen = "NaN"
	return sen

def tt_to_par_tt_conversion(par, value):
	# Converts the tt as an integer, to the tt measurement made in 
	# the PAR  data file 
	temp_par_df = par[par.tt > value]
	temp_tt = temp_par_df.iloc[0][0]
	return temp_tt

def circ_days(par, start, stop):
	# Calculates the number of days between two tt index values
	temp_start = par[par.tt == start]
	temp_stop = par[par.tt == stop]
	return temp_stop - temp_start

def daily_rfr_calc(par_index, par, rfr_list):
	# Find the PAR interception on a given day
	tt = par.iloc[par_index][0]
	tt = int(round(tt))
	daily_rfr = rfr_list[(a + 1)]
	return daily_rfr

def daily_par_calc(par, daily_rfr, par_index):
	# Calculate the PAR intercepted on a given day
	daily_par = daily_rfr * (0.5 * par["kipp"][par_index])	
	return daily_par

def total_par_calc(start, stop):
	# Calculate the PAR intercepted between two days

	# Convert the input termal time values, into their nearest 
	# equivalents in the par file. The corresponding
	# par measurements can then be used 
	par_tt_start = tt_to_par_tt_conversion(par, start)
	par_tt_stop = tt_to_par_tt_conversion(par, stop)
	
	canopy_duration_circ = circ_days(par, par_tt_start, par_tt_stop)
	
	daily_par_list = []
	for count in range(0, canopy_duration_circ):
		daily_rfr = daily_rfr_calc((par_index + count))
		daily_par = daily_par_calc(par, daily_rfr, (par_index + count))
		daily_par_list.append(daily_par)

	total_par = sum(daily_par_list)
	return total_par, circ_days
	
df = data_input("C:\\users\\john\\google drive\\modelling\\raw.csv")
