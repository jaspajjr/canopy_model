import time
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import sympy as sym 
from sympy import exp
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
		sen = np.NAN
		return sen
	else:
		sen = temp_sen_list[-1]
		if sen < 1000:
			sen = np.NAN 
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

def par_calc(start, stop, par, plot, p0):
	'''calculate the time par itercepted between a and b'''
	# Need to get tt value from every day between start and stop
	# Calculate the index of start
	if np.isnan(start) == True or np.isnan(stop) == True:
		return np.NAN, np.NAN
	start_index = par[par.tt > start].index[0]
	# Calculate the index of stop
	stop_index = par[par.tt > stop].index[0]
	# Calculate duration in circadian days
	circ_days = stop_index - start_index
	# Need to get p0 and tt 
	p0_plot = p0.loc["Plot%d" %plot]
	# Need to get rfr value for each tt between start and stop
	x = [par["tt"].loc[count] for count in range(start_index, stop_index)]
	# Need to get KIPP value from every day between start and stop
	kipp = [par["kipp"].loc[count] for count in range(start_index, stop_index)]
	# Calculate f(x) for each x in the x list
	rfr = [f(item, p0_plot) for item in x]
	# Calculate PAR value from the kipp values
	par = [(0.5 * item) for item in kipp]
	# Convert list to Series
	rfr = pd.Series(rfr)
	# Multiply each rfr value by corresponding KIPP / 2
	par = pd.Series(par)
	par_int = rfr * par
	# Sum the results of that for each tt 
	t_par = sum(par_int)
	return t_par, circ_days

def rfr_n(start, n, plot):
	''' Calculate the change in rfr value n days after 
	reference '''
	if np.isnan(start) == True:
		return np.NAN
	# Find the index of the starting tt value
	start_index = par[par.tt > start].index[0]
	# get the index n days after start_index
	out_index = start_index + n
	out_tt = par["tt"].ix[out_index]
	p0_plot = p0.loc["Plot%d" %plot]
	rfr_out = f(out_tt, p0_plot)
	return rfr_out

'''
Derivative related modules
'''
def f_prime(reference_x, n, plot):
	''' Calculate the rate of change in the canopy formation at 
	point reference_x, given the parameters for the given plot in p0'''
	if np.isnan(reference_x) == True:
		return np.NAN
	# Find index of the starting tt value
	start_index = par[par.tt > reference_x].index[0]
	# get the index n days after start_index
	out_index = start_index + n
	out_tt = par["tt"].ix[out_index]
	if plot == 5:
		print out_tt
		print

	p0_plot = p0.loc["Plot%d" %plot]
	if np.isnan(p0_plot[0]) == True:
		return np.NAN
	else:
		x = sym.symbols("x")
		c, b1, m1, a, b2, m2 = p0_plot
		f = (c / (1 + exp(-b1 * (x - m1)))) * (a + (c * exp(-exp(-b2 * (x - m2)))))
		fx = sym.diff(f, x)
		dx_dy = fx.evalf(n=35, subs={x: out_tt})
	return dx_dy
'''
Integral related modules
'''
def int_calc(start, stop, plot):
	# Calculates the integral of the function f_temp, the main function, with the
	# parameters of the specific plot, between points start and stop
	if np.isnan(start) == True or np.isnan(stop) == True:
		return np.NAN, np.NAN
	else:
		p0_plot = p0.loc["Plot%d" %plot]
		c, b1, m1, a, b2, m2 = p0_plot 
		f_temp = lambda x: (c / (1 + np.exp(-b1 * (x - m1)))) * (a + (c * np.exp((-1 * (np.exp(-b2 * (x - m2)))))))
		ans, err = quad(f_temp, start, stop)
	return ans, err 

'''
########################################################################
This section of the code deals with functions being executed, rather than
functions being created. 
'''
''''''''''''''''''''''''''''''''''''''''''''''''
# Input the raw data
#df = pd.read_csv("C:\\users\\john\\google drive\\modelling\\raw.csv")
df = pd.read_csv("C:\\users\\john\\google drive\\modelling\\raw_allele.csv")
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
	p0["Plot%d" %(plot + 1)] = temp_p0
marquardt_time = time.time() - marquardt_start
p0 = p0.transpose()
p0.columns = ["c", "b1", "m1", "a", "b2", "m2"]
print "Marquardt time %d" %marquardt_time
#p0.to_csv("C:\\users\\john\\canopy_model\\p0.csv")
p0.to_csv("C:\\users\\john\\canopy_model\\p0_allele.csv")

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
#rfr_df.to_csv("C:\\users\\john\\canopy_model\\rfr.csv")
rfr_df.to_csv("C:\\users\\john\\canopy_model\\rfr_allele.csv")
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
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate PAR between sowing and senescence
t_par_start = time.time()
t_par = []
t_par_dur = []
plot_count = 1
for item in field["sen"]:
	if np.isnan == True:
		temp_par_val = np.NAN, np.NAN
	else:
		temp_par_val = par_calc(0, item, par, plot_count, p0)
	t_par.append(temp_par_val[0])
	t_par_dur.append(temp_par_val[1])
	plot_count += 1
field["t_par"] = t_par
field["canopy_duration"] = t_par_dur
t_par_elapsed = time.time() - t_par_start
print "Total PAR intercepted %d seconds" %t_par_elapsed
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the PAR between sowing and GS31
gs31_par_time = time.time()
gs31_par = []
gs31_par_dur = []
plot_count = 1
for stop in field["gs31"]:
	if np.isnan(stop) == True:
		temp_par_val = np.NAN, np.NAN
	else:
		temp_par_val = par_calc(0, stop, par, plot_count, p0)
	gs31_par.append(temp_par_val[0])
	gs31_par_dur.append(temp_par_val[1])
	plot_count += 1
field["gs31_par"] = gs31_par
field["sow_gs31_dur"] = gs31_par_dur
gs31_par_elapsed = time.time() - gs31_par_time
print "PAR to GS31 %d" %gs31_par_elapsed
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate PAR between gs31_anth
gs31_anth_par_time = time.time()
gs31_anth_par_list = []
gs31_anth_par_dur = []
plot_count = 1
for start, stop in zip(field["gs31"], field["anth"]):
	if np.isnan(start) == True or np.isnan(stop) == True:
		temp_par_val = np.NAN, np.NAN
	else:
		temp_par_val = par_calc(start, stop, par, plot_count, p0)
	gs31_anth_par_list.append(temp_par_val[0])
	gs31_anth_par_dur.append(temp_par_val[1])
	plot_count += 1
field["par_gs31_anth"] = gs31_anth_par_list
field["gs31_anth_dur"] = gs31_anth_par_dur
par_gs31_anth_elapsed = time.time() - gs31_anth_par_time
print "GS31 to Anth PAR time %d" %par_gs31_anth_elapsed

''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the par between anth and sen
anth_sen_par_time = time.time()
anth_sen_par_list = []
anth_sen_par_dur_list = []
plot_count = 1
for start, stop in zip(field["anth"], field["sen"]):
	if np.isnan(start) == True or np.isnan(stop) == True:
		temp_par_val = np.NAN, np.NAN
	else:
		temp_par_val = par_calc(start, stop, par, plot_count, p0)
	anth_sen_par_list.append(temp_par_val[0])
	anth_sen_par_dur_list.append(temp_par_val[1])
	plot_count += 1
field["par_anth_sen"] = anth_sen_par_list
field["anth_sen_dur"] = anth_sen_par_dur_list
par_anth_sen_elapsed = time.time() - anth_sen_par_time
print "Anth - Sen PAR time %d seconds" %par_anth_sen_elapsed

''''''''''''''''''''''''''''''''''''''''''''''''
rfr_n_time = time.time()
# Calculate the rfr values 5, 10, 15 days after anthesis
field["rfr_5_anth"] = [rfr_n(field["anth"].iloc[(x - 1)], 5, x) for x in xrange(1, (len(field) + 1))]
field["rfr_10_anth"] = [rfr_n(field["anth"].iloc[(x - 1)], 10, x) for x in xrange(1, (len(field) + 1))]
field["rfr_15_anth"] = [rfr_n(field["anth"].iloc[(x - 1)], 15, x) for x in xrange(1, (len(field) + 1))]
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the rfr value 5, 10, 15 days before anthesis
field["rfr_5_before_anth"] = [rfr_n(field["anth"].iloc[(x - 1)], -5, x) for x in xrange(1, (len(field) + 1))]
field["rfr_10_before_anth"] = [rfr_n(field["anth"].iloc[(x - 1)], -10, x) for x in xrange(1, (len(field) + 1))]
field["rfr_15_before_anth"] = [rfr_n(field["anth"].iloc[(x - 1)], -15, x) for x in xrange(1, (len(field) + 1))]
rfr_n_elapsed = time.time() - rfr_n_time
print "R:FR after n completed, %d" %rfr_n_elapsed
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the integral between the start and end of the canopy

int_sen_time = time.time()
int_list = []
plot_count = 1
for stop in field["sen"]:
	if np.isnan(stop) == True:
		temp_int_val = np.NAN, np.NAN
	else:
		temp_int_val = int_calc(0, stop, plot_count)
	int_list.append(temp_int_val[0])
	plot_count += 1
field["int_sen"] = int_list
int_sen_elapsed = time.time() - int_sen_time
print "Integral from start to senescence %d" %int_sen_elapsed
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the integral between the start and gs31

int_gs31_time = time.time()
int_gs31_list = []
plot_count = 1
for stop in field["gs31"]:
	if np.isnan(stop) == True:
		temp_int_val = np.NAN, np.NAN
	else:
		temp_int_val = int_calc(0, stop, plot_count)
	int_gs31_list.append(temp_int_val[0])
	plot_count += 1
field["int_gs31"] = int_gs31_list
int_gs31_elapsed = time.time() - int_gs31_time
print "integral to gs31 %d" %int_gs31_elapsed
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the integral between gs31 and anth

int_gs31_anth_time = time.time()
int_gs31_anth_list = []
plot_count = 1
for start, stop in zip(field["gs31"], field["anth"]):
	if np.isnan(start) == True or np.isnan(stop) == True:
		temp_int_val = np.NAN, np.NAN
	else:
		temp_int_val = int_calc(start, stop, plot_count)
	int_gs31_anth_list.append(temp_int_val[0])
	plot_count += 1
field["int_gs31_anth"] = int_gs31_anth_list
int_gs31_anth_elapsed = time.time() - int_gs31_anth_time
print "integral between gs31 and anth %d" %int_gs31_anth_elapsed
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the integral between anthesis and senescence
int_anth_sen_time = time.time()
int_anth_sen_list = []
plot_count = 1
for start, stop in zip(field["anth"], field["sen"]):
	if np.isnan(start) == True or np.isnan(stop) == True:
		temp_int_val = np.NAN, np.NAN
	else:
		temp_int_val = int_calc(start, stop, plot_count)
	int_anth_sen_list.append(temp_int_val[0])
	plot_count += 1
field["int_anth_sen"] = int_anth_sen_list
int_anth_sen_elapsed = time.time() - int_anth_sen_time
print "integral between anthesis and senescence %d" %int_anth_sen_elapsed
''''''''''''''''''''''''''''''''''''''''''''''''
dy_time = time.time()
# Calculate the derivative at gs31 and at anth
field["dy_at_gs31"] = [f_prime(field["gs31"].iloc[(x - 1)], 0, x) for x in xrange(1, (len(field) + 1))]
field["dy_at_anth"] = [f_prime(field["anth"].iloc[(x - 1)], 0, x) for x in xrange(1, (len(field) + 1))]
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the derivative 5, 10, 15 days after anthesis
field["dy_5_post_anth"] = [f_prime(field["anth"].iloc[(x - 1)], 5, x) for x in xrange(1, (len(field) + 1))]
field["dy_10_post_anth"] = [f_prime(field["anth"].iloc[(x - 1)], 10, x) for x in xrange(1, (len(field) + 1))]
field["dy_15_post_anth"] = [f_prime(field["anth"].iloc[(x - 1)], 15, x) for x in xrange(1, (len(field) + 1))]
''''''''''''''''''''''''''''''''''''''''''''''''
# Calculate the derivative 5, 10, 15 days before anthesis
field["dy_5_pre_anth"] = [f_prime(field["anth"].iloc[(x - 1)], -5, x) for x in xrange(1, (len(field) + 1))]
field["dy_10_pre_anth"] = [f_prime(field["anth"].iloc[(x - 1)], -10, x) for x in xrange(1, (len(field) + 1))]
field["dy_15_pre_anth"] = [f_prime(field["anth"].iloc[(x - 1)], -15, x) for x in xrange(1, (len(field) + 1))]


dy_time_elapsed = time.time() - dy_time
print "Derivative time %d" %dy_time_elapsed
''''''''''''''''''''''''''''''''''''''''''''''''

#field.to_csv("C:\\users\\john\\google drive\\modelling\\canopy_model.csv")
field.to_csv("C:\\users\\john\\google drive\\modelling\\canopy_model_allele.csv")

total_time = time.time() - start_time
print "Modelling finished"
print "Total time taken %d" %total_time
