import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

# Input data

df = pd.read_csv("C:\\users\\john\\canopy_model\\rfr.csv")
line_names = pd.read_csv("C:\\users\\john\\canopy_model\\line_names.csv")

df["Line"] = line_names["Line"]
# ID'd lines
ID = ["Ren", "Sav", "SRen103", "SRen104", "SRen4", "SRen51", "SRen76", "SRen61"]

# Mean of genotype
df = df.groupby("Line").mean()


''''''''''''''''''''''''''''''
# Plotting
def scatter(ID, df):
	x = [int(item) for item in df.columns]
	color_list = ["blue", "green", "red", "cyan", "magenta", "yellow", "black", "burlywood"]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	axis_text_size = 8
	for item in range(0, len(ID)):
		ax.plot(x, df.ix[ID[item]], color=color_list[item], linestyle='-', label=ID[item])
	ax.set_xlabel("Thermal time ($degree$ $C$)")
	ax.set_ylabel("R:FR")
	ax.legend(loc=8, prop={"size": axis_text_size})
	ax.tick_params(axis='both', which='major', labelsize=axis_text_size)
	ax.xaxis.grid(color='gray', linestyle="--")
	ax.yaxis.grid(color='gray', linestyle='--')
	ax.set_xlim(0, 3000)
	ax.set_ylim(0, 1)
	plt.show()

scatter(ID, df)

'''
fig = plt.figure(111)
def scatter(x_series, y_series, descriptions_dict):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	axis_text_size = 8
	#x_series = x_series.dropna()
	#y_series = y_series.dropna()
	x_list = x_series.loc["Ren":"SRen99"]
	y_list = y_series.loc["Ren":"SRen99"]
	x_err = x_series.loc["Avg"]
	y_err = y_series.loc["Avg"]
	ax.plot(x_list, y_list, color='black', marker='x', linestyle='none',
		markersize=8)
	ax.set_xlabel(descriptions_dict[(x_series.name)])
	ax.set_ylabel(descriptions_dict[y_series.name])
	#ax.set_title(r'%s $\alpha \degree > \beta $' %test)
	ax.xaxis.grid(color='gray', linestyle='--')
	ax.yaxis.grid(color="gray", linestyle='--')
	ax.tick_params(axis='both', which='major', labelsize=axis_text_size)
	annotate_list = annotate_maker(x_series, y_series) 
	for text, x, y in annotate_list:
		ax.annotate(text, xy=(x, y), 
			xytext=(-40, 20), textcoords='offset points', 
			arrowprops=dict(arrowstyle="->", connectionstyle='arc3, rad=0.5',
			color='gray'), fontsize=8)
	#ax.autoscale(enable=True, tight=True, axis=x)
'''