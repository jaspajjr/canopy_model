import time
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import sympy as sym 
from scipy.optimize import curve_fit
from scipy.integrate import quad
start_time = time.time()
print "Starting"

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
		out = [np.NAN, np.NAN, np.NAN, np.NAN, np.NAN, np.NAN]
	return out

def plot_rfr(p0):
	if type(p0[0]) == str:
		y_pred = np.NAN
		return y_pred
	y_pred = [f(a, p0) for a in xrange(0, 3001)]
	return y_pred

def max_rfr_finder(rfr_list):
	# Finds the maximum r:fr and the tt to the max r:fr value 
	max_rfr = max(rfr_list)
	max_rfr_tt = rfr_list.index(max_rfr)
	return max_rfr, max_rfr_tt

def senescence(max_rfr, sen_rfr, rfr_list):
	'''Perform type check before this'''
	# Calculates the tt to senescence
	rfr_list = rfr_list.tolist()
	#for item in (rfr_list[(rfr_list.index(round(max_rfr, 2))):]):
	count = 0
	temp_sen_list = []
	for item in rfr_list:
		count += 1
		if round(item, 2) == round(sen_rfr, 2):
			temp_sen_list.append(count)

	if len(temp_sen_list) == 0:
		sen = "NaN"
		return sen
	else:
		sen = temp_sen_list[-1]
	return sen

def tt_to_par_tt_conversion(par, value):
	# Converts the tt as an integer, to the tt measurement made in 
	# the PAR  data file, outputs as float64 to allow for recognition
	# in the PAR calculation modules 
	try:
		temp_par_df = par[par.tt > value]
		temp_tt = np.float64(temp_par_df.iloc[0][0])
	except IndexError:
		temp_tt = np.float64(0)
	return temp_tt

def tt_to_par_tt_rfr_n(par, value):
	# Converts the tt as an integer, to the tt measurement made in 
	# the PAR  data file, outputs as float64 to allow for recognition
	# in the PAR calculation modules 
	try:
		temp_par_df = par[par.tt > value]
		temp_tt = int(temp_par_df.iloc[0][0])
	except IndexError:
		temp_tt = int(0)
	return temp_tt

def circ_days(par, start, stop):
	# Calculates the number of days between two tt index values
	if start == 0:
		temp_start = 0
	else:
		try:
			temp_start = par[par.tt == start].index[0]
		except IndexError:
			print "Shit, index Error"
			print start
	if stop == 0:
		return np.NAN
	else:
		temp_stop = par[par.tt == stop].index[0]
	return temp_stop - temp_start

def daily_rfr_calc(par_index, par, rfr_list):
	# Find the PAR interception on a given day
	tt = par.iloc[par_index][0]
	tt = int(round(tt))
	daily_rfr = rfr_list[(par_index + 1)]
	return daily_rfr

def daily_par_calc(par, daily_rfr, par_index):
	# Calculate the PAR intercepted on a given day
	daily_par = daily_rfr * (0.5 * par["kipp"][par_index])	
	return daily_par

def total_par_calc(start, stop, par, plot):
	# Calculate the PAR intercepted between two days

	# Convert the input termal time values, into their nearest 
	# equivalents in the par file. The corresponding
	# par measurements can then be used 
	par_tt_start = tt_to_par_tt_conversion(par, start)
	par_tt_stop = tt_to_par_tt_conversion(par, stop)
	if np.isnan(par_tt_stop) == True:
		total_par = np.NAN
		circ_days_duration = np.NAN
		return total_par, circ_days_duration	
	else:
		canopy_duration_circ = circ_days(par, par_tt_start, par_tt_stop)
		daily_par_list = []
		for par_index in range(par_tt_start, par_tt_stop):
			try:
				daily_rfr = daily_rfr_calc((par_index), par, rfr_df.iloc[plot])
				daily_par = daily_par_calc(par, daily_rfr, (par_index))
				daily_par_list.append(daily_par)
			except IndexError:
				daily_rfr = 0
		
		total_par = sum(daily_par_list)	
	return total_par, canopy_duration_circ

def rfr_after_n(reference, n_days, par, plot):
	# Calculates the absolute R:FR value n days after the reference
	# point. n_days is expected to be an integer
	if np.isnan(reference) == True:
		rfr_after_n_days = np.NAN
	else:
		reference_index = tt_to_par_tt_rfr_n(par, reference)
		n_days = reference_index + n_days
		rfr_list = rfr_df.loc["Plot%d" %plot].tolist()
		n_days = par[par.tt > n_days].iloc[n_days][0]
		rfr_after_n_days = daily_rfr_calc(n_days, par, rfr_list)
	return rfr_after_n_days

'''
Derivative related modules
'''
'''def dy_f(p0, x):
	# The derivative of the function f(x) at x
	c, b1, m1, a, b2, m2 = p0
	dy_f_x = ((b1 * c * exp(-b1 * (-m1 + x))) / (np.power(1 + exp(-b1 * (-m1 + x)), 2))) * (c * exp(-exp(-b2 * (-m2 + x))) + a) + ((np.power(c, 2) * b2 * exp(-b2 * (-m2 + x)) * exp(-exp(-b2 * (-m2 + x)))) / 1 + exp(-b1 * (-m1 + x)))
	return dy_f_x

'''

'''def f_prime(x, p0):
	#Calculates the value of the derivative according to p0
	at the given x value
	
	c, b1, m1, a, b2, m2 = p0'''

'''
Integral related modules
'''

def integral_calc(p0, start, stop):
	''' Calculates the definite integral between start and stop, 
	given parameters p0
	'''
	# The equation to be fitted, however with the parameters stored in p0
	c, b1, m1, a, b2, m2 = p0
	temp_f = (c / (1 + np.exp(-b1 * (x - m1)))) * (a + (c * np.exp((-1 * (np.exp(-b2 * (x - m2)))))))
	ans, err = quad(temp_f, start, stop)
	return ans

'''
########################################################################
This section of the code deals with functions being executed, rather than
functions being created. 
'''
''''''''''''''''''''''''''''''''''''''''''''''''
# Input the raw data
df = pd.read_csv("C:\\users\\john\\google drive\\modelling\\raw.csv")
par = pd.read_csv("C:\\users\\john\\google drive\\modelling\\par.csv")
anth = pd.read_csv("C:\\users\\john\\google drive\\modelling\\anth.csv")

tt = df.loc[0]
df = df.loc[1:]
''''''''''''''''''''''''''''''''''''''''''''''''
# Set up the dataframe which will ouput the final values
field = pd.DataFrame(index=df.index)
field["gs31"] = anth["gs_31"]
field["anth"] = anth["anth"]
field["gs31_anth"] = field["anth"] - field["gs31"]

''''''''''''''''''''''''''''''''''''''''''''''''
# Setup the dataframe which will store the parameter values
p0_index = [x for x in xrange(0, 6)]
p0 = pd.DataFrame(index=p0_index)
# Remove the anth dataframe as it is no longer needed
del anth

''''''''''''''''''''''''''''''''''''''''''''''''
# Setup for the loop to calculate p0 from marquardt
p0_initial = [0.8, 0.001, 2000, 0, -0.001, 700]
marquardt_start = time.time()
for plot in xrange(0, len(df)):
	temp = []
	temp_p0 = marquardt(f_fit, tt, df.iloc[plot], p0_initial)
	p0["Plot%d" %plot] = temp_p0
marquardt_time = time.time() - marquardt_start
p0 = p0.transpose()
p0.columns = ["c", "b1", "m1", "a", "b2", "m2"]
print "Marquardt time %d" %marquardt_time

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculating the residuals
res_index = [x for x in xrange(0, 21)]
residuals = pd.DataFrame(index=res_index)
res_start = time.time()
for plot in xrange(1, len(df)):
	y_est, res2 = residual_calculator(tt, df.loc[1], p0.iloc[(plot-1)])
	residuals["Plot%d" %plot] = res2
res_time = time.time() - res_start
residuals = residuals.transpose()
cols = ["r%d" %x for x in xrange(1, 22)]
print "Residual time %d" %res_time

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate rfr values for each plot
rfr_df_index = [x for x in xrange(0, 3001)]
rfr_df = pd.DataFrame(index=rfr_df_index)
rfr_start = time.time()
#rfr_df["Plot0"] = [0 for x in xrange(0, 3001)]
for plot in xrange(0, (len(df))):
	rfr_df["Plot%d" %(plot + 1)] = plot_rfr(p0.iloc[(plot)])
rfr_df = rfr_df.transpose()
rfr_time = time.time() - rfr_start
print "rfr time %d" %rfr_time

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the maximum rfr, and tt to max rfr
max_rfr_finder_time_start = time.time()
max_rfr_list = []
max_rfr_tt_list = []
for plot in xrange(0, len(df)):
	if type(rfr_df.iloc[plot][0]) == str:
		max_rfr_temp, max_rfr_tt_temp = np.NAN, np.NAN
	else:
		max_rfr_temp, max_rfr_tt_temp = max_rfr_finder(rfr_df.iloc[plot].tolist())
	max_rfr_list.append(max_rfr_temp)
	max_rfr_tt_list.append(max_rfr_tt_temp)
field["max_rfr"] = max_rfr_list
field["max_rfr_tt"] = max_rfr_tt_list
max_rfr_finder_time = time.time() - max_rfr_finder_time_start
print "Max R:FR time %d" %max_rfr_finder_time

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculates the 37% of the maximum R:FR
field = field.convert_objects(convert_numeric=True)
field["sen_mult"] = [0.37 for x in xrange(0, len(df))]
field["sen_rfr"] = field["max_rfr"] * field["sen_mult"]
del field["sen_mult"]

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate time to senescence
sen = []
sen_start = time.time()
count = 0
for plot in xrange(0, (len(df))):
	count += 1
	if type(rfr_df.loc["Plot%d" %(plot + 1)][0]) == str:
		sen_temp = np.NAN
	else:
		sen_temp = senescence((field.iloc[plot]["max_rfr"]), 
			field.iloc[plot]["sen_rfr"],
			rfr_df.loc["Plot%d" %(plot + 1)])
	sen.append(sen_temp)
field["sen"] = sen
field = field.convert_objects(convert_numeric=True)
senescence_time = time.time() - sen_start
print "Senescence time %d" % senescence_time

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate time from anthesis to senescence

field["anth_sen"] = field["sen"] - field["anth"]

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the total PAR intercepted over the entire season
# This could have been done as an apply function using the pandas
# library, however, at the time of writing I didn't plan for that 
# so it is quicker and more efficient to use the code as intended 
# than to rewrite it and make it more pythonic
t_par_start = time.time()
temp_par_list = []
t_par_can_dur = []
plot_count = 1
for item in field["sen"]:
	if np.isnan(item) == True:
		temp_par_val = np.NAN, np.NAN
	else:
		temp_par_val = total_par_calc(0, item, par, plot_count)
	print temp_par_val[0]
	temp_par_list.append(temp_par_val[0])
	t_par_can_dur.append(temp_par_val[1])
	plot_count += 1
print sum(temp_par_val)	
field["total_par"] = temp_par_list
field["t_can_dur_circ"] = t_par_can_dur
t_par_elapsed = time.time() - t_par_start
print "Total PAR time %d" %t_par_elapsed

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the PAR intercepted between sowing and gs31
gs31_par_start = time.time()
temp_gs31_par_list = []
temp_gs31_par_dur = []
plot_count = 1
for item in field["gs31"]:
	if np.isnan(item) == True:
		temp_par_val = np.NAN, np.NAN
	else:
		temp_par_val = total_par_calc(0, item, par, plot_count)
	temp_gs31_par_list.append(temp_par_val[0])
	temp_gs31_par_dur.append(temp_par_val[1])
	plot_count += 1	
field["par_gs31"] = temp_gs31_par_list
field["par_gs31_dur_circ"] = temp_gs31_par_dur
par_gs31_elapsed = time.time() - gs31_par_start
print "GS31 PAR time %d" %par_gs31_elapsed

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate PAR between gs31_anth
gs31_anth_par_start = time.time()
temp_gs31_anth_par_list = []
temp_gs31_anth_par_dur = []
plot_count = 1
for start, stop in zip(field["gs31"], field["anth"]):
	if np.isnan(start) == True or np.isnan(stop) == True:
		temp_par_val = np.NAN, np.NAN
	else:
		temp_par_val = total_par_calc(start, stop, par, plot_count)
	temp_gs31_anth_par_list.append(temp_par_val[0])
	temp_gs31_anth_par_dur.append(temp_par_val[1])
	plot_count += 1
field["par_gs31_anth"] = temp_gs31_anth_par_list
field["par_gs31_anth_dur_circ"] = temp_gs31_anth_par_dur
par_gs31_anth_elapsed = time.time() - gs31_anth_par_start
print "GS31 to Anth PAR time %d" %par_gs31_anth_elapsed

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the par between anth and sen
anth_sen_par_time = time.time()
temp_anth_sen_par_list = []
temp_anth_sen_par_dur = []
plot_count = 1
for start, stop in zip(field["anth"], field["sen"]):
	if np.isnan(start) == True or np.isnan(stop) == True:
		temp_par_val = np.NAN, np.NAN
	else:
		temp_par_val = total_par_calc(start, stop, par, plot_count)
	temp_anth_sen_par_list.append(temp_par_val[0])
	temp_anth_sen_par_dur.append(temp_par_val[1])
	plot_count += 1
field["par_anth_sen"] = temp_anth_sen_par_list
field["par_anth_sen_dur_circ"] = temp_anth_sen_par_dur
par_anth_sen_elapsed = time.time() - anth_sen_par_time
print "PAR Anth - Sen %d" %par_anth_sen_elapsed

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the rfr 5 days after anthesis
#field["rfr_5_after_anth"] = 

field.to_csv("C:\\users\\john\\google drive\\modelling\\canopy_model_test.csv")

print field.head()
total_time = time.time() - start_time
print "Total time taken %d" %total_time
a = [rfr_after_n(field["anth"][x], 5, par, (x)) for x in range(1, (len(field) + 1))]
