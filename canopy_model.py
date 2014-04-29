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