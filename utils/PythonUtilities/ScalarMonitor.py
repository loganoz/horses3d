#!/usr/bin/env python
#####################
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def GatherNewScalarValues(fileName,skip_data):
	fid = open(fileName, 'r')
	counter = -1
	iter = []
	time = []
	value = []
#
#	First line is the monitor name
#	------------------------------
	line = fid.readline()
	entries = line.split()
	monitorName = entries[-1]
#
#	Navigate until the data beginning
#	---------------------------------
	for line in fid:
		entries = line.split()
		if ( len(entries) == 0 ):
			continue
		if ( "Selected variable:" in line ):
			variable = line.replace("Selected variable: ","");
			variable = variable.lower()
		elif ( entries[0] == "Iteration" ):
			break 
#
#	Now skip the requested amount of data
#	-------------------------------------
	if ( skip_data > 0 ):
		for i in xrange(int(skip_data)-1):
			fid.next()

	for line in fid:
		numbers = line.split()
		iter.append(float(numbers[0]))
		time.append(float(numbers[1]))
		value.append(float(numbers[2]))
		counter = counter + 1

	if ( counter == -1 ):
		counter = 0

	return (counter,time,value, monitorName, variable)

def NewScalarPlot(time, value, ax, monitorName,variable):
	if ( "kinetic energy rate" in variable ):
		ax.plot(time,[x * (-1) for x in value],'-',color='#1f77b4', label=monitorName, linewidth=1.5)
	else:
		ax.plot(time, value,'-',color='#1f77b4', label=monitorName, linewidth=1.5)

	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles,labels)

def AppendToScalarPlot(time, value, ax,variable):
	if ( "kinetic energy rate" in variable ):
		ax.plot(time,[x * (-1) for x in value],'-',color='#1f77b4',linewidth=1.5)
	else:
		ax.plot(time, value,'-',color='#1f77b4',linewidth=1.5)


def UpdateScalarPlot(fileName, skip_data, ax):
	N,time,value,monitorName,variable = GatherNewScalarValues(fileName, skip_data)

	if ( skip_data == 0 ):
		NewScalarPlot(time, value, ax, monitorName, variable)
	else:
		AppendToScalarPlot(time, value, ax, variable)


	return (N)
