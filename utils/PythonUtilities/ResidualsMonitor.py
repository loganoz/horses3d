#!/usr/bin/env python
#####################
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

RESIDUALS_HEADER_SIZE = 2

def GatherNewResidualsValues(fileName,skip_data):
	residuals_fid = open(fileName, 'r')
	counter = -1
	iter = []
	time = []
	continuity = []
	x_momentum = []
	y_momentum = []
	z_momentum = []
	energy = []

	for i in range(RESIDUALS_HEADER_SIZE):
		next(residuals_fid)

	if ( skip_data > 0 ):
		for i in range(int(skip_data)-1):
			next(residuals_fid)

	for line in residuals_fid:
		numbers = line.split()
		iter.append(float(numbers[0]))
		time.append(float(numbers[1]))
		continuity.append(float(numbers[3]))
		x_momentum.append(float(numbers[4]))
		y_momentum.append(float(numbers[5]))
		z_momentum.append(float(numbers[6]))
		energy.append(float(numbers[7]))
		counter = counter + 1

	if ( counter == -1 ):
		counter = 0

	return (counter,time,continuity,x_momentum,y_momentum,z_momentum,energy)

def NewResidualsPlot(time, continuity, x_momentum, y_momentum, z_momentum, energy, ax):
	ax.semilogy(time,continuity,'-',color='#1f77b4', label='continuity',linewidth=1.5)
	ax.semilogy(time,x_momentum,'-',color='#ff7f0e', label='x-momentum',linewidth=1.5)
	ax.semilogy(time,y_momentum,'-',color='#2ca02c', label='y-momentum',linewidth=1.5)
	ax.semilogy(time,z_momentum,'-',color='#d62728', label='z-momentum',linewidth=1.5)
	ax.semilogy(time,energy,'-',color='#bcbd22', label='energy',linewidth=1.5)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles,labels)
	ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=5, mode="expand", borderaxespad=0.,prop={'size': 10})

def AppendToResidualsPlot(time, continuity, x_momentum, y_momentum, z_momentum, energy, ax):
	ax.semilogy(time,continuity,'-',color='#1f77b4',linewidth=1.5)
	ax.semilogy(time,x_momentum,'-',color='#ff7f0e',linewidth=1.5)
	ax.semilogy(time,y_momentum,'-',color='#2ca02c',linewidth=1.5)
	ax.semilogy(time,z_momentum,'-',color='#d62728',linewidth=1.5)
	ax.semilogy(time,energy,'-',color='#bcbd22',linewidth=1.5)


def UpdateResidualsPlot(fileName, skip_data, ax):
	N,time,continuity,x_momentum,y_momentum,z_momentum,energy = GatherNewResidualsValues(fileName, skip_data)

	if ( skip_data == 0 ):
		NewResidualsPlot(time, continuity, x_momentum, y_momentum, z_momentum, energy, ax)
	else:
		AppendToResidualsPlot(time, continuity, x_momentum, y_momentum, z_momentum, energy, ax)


	return (N)
