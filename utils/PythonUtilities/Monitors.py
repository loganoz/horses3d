#!/usr/bin/env python
#####################
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from ResidualsMonitor import *
from ScalarMonitor import *

RESIDUALS_FILE = 1
MONITOR_FILE = 2

def getWhichMonitorType(fileName):
	if ".residuals" in fileName:
		typeIs = RESIDUALS_FILE
	else:
		typeIs = MONITOR_FILE

	return (typeIs)
	
def plotMonitors(fileNames):

	if ( isinstance(fileNames,str) ):
		fileNames = [fileNames]

	no_of_plots = len(fileNames)
	skip_data = np.zeros(no_of_plots)
	plotType = [];
		
	for i in range(0, no_of_plots):
		plotType.append(getWhichMonitorType(fileNames[i]))
#	
#	Initialize plt
#	--------------
	h = plt.figure(facecolor="white")

	plt.ion()

	while True:
		for i in range(0,no_of_plots):
			ax = plt.subplot(no_of_plots,1,i+1)
			ax.set_facecolor((0.05,0.05,0.05))
			plt.grid(b=True, which='both', color='0.2',linestyle='-')

			if ( plotType[i] == RESIDUALS_FILE ):
				N = UpdateResidualsPlot(fileNames[i], skip_data[i], ax)		
			elif ( plotType[i] == MONITOR_FILE ):
				N = UpdateScalarPlot(fileNames[i], skip_data[i], ax)		

			skip_data[i] = skip_data[i] + N

		plt.pause(2)
