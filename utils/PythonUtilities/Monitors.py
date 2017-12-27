#!/usr/bin/env python
#####################
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from ResidualsMonitor import *

RESIDUALS_FILE = 1
MONITOR_FILE = 2

fileName = ["/Users/juanmanzanero/OwnCloud/Research/DGSEM/Codes/FORTRAN/HORSES3D/Solver/test/NavierStokes/TaylorGreen/RESULTS/TGV.residuals"]

def getWhichMonitorType(fileName):
	if ".residuals" in fileName:
		typeIs = RESIDUALS_FILE
	else:
		typeIs = MONITOR_FILE

	return (typeIs)
	
def plotMonitors(fileNames):

	no_of_plots = len(fileNames)
	skip_data = np.zeros(no_of_plots)
	plotType = [];
		
	for i in range(0, no_of_plots):
		plotType.append(getWhichMonitorType(fileNames[i]))
	

	h = plt.figure(1);

	plt.ion()

	while True:
		for i in range(0,no_of_plots):
			if ( plotType[i] == RESIDUALS_FILE ):
				ax = plt.subplot(no_of_plots,1,i+1)
				N = UpdateResidualsPlot(fileNames[i], skip_data[i], ax)		
				skip_data[i] = skip_data[i] + N

		plt.pause(2)


plotMonitors(fileName)

files = ["a","b","c"];
print(len(files))
