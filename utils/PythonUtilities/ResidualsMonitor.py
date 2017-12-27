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
	fig = plt.figure(1)
	counter = -1
	iter = []
	time = []
	continuity = []
	x_momentum = []
	y_momentum = []
	z_momentum = []
	energy = []

	for i in xrange(RESIDUALS_HEADER_SIZE):
		residuals_fid.next()

	if ( skip_data > 0 ):
		for i in xrange(int(skip_data)-1):
			residuals_fid.next()

	for line in residuals_fid:
		numbers = line.split()
		iter.append(float(numbers[0]))
		time.append(float(numbers[1]))
		continuity.append(float(numbers[2]))
		x_momentum.append(float(numbers[3]))
		y_momentum.append(float(numbers[4]))
		z_momentum.append(float(numbers[5]))
		energy.append(float(numbers[6]))
		counter = counter + 1

	if ( counter == -1 ):
		counter = 0

	return (counter,time,continuity,x_momentum,y_momentum,z_momentum,energy)

def NewResidualsPlot(time, continuity, x_momentum, y_momentum, z_momentum, energy, ax):
	ax.semilogy(time,continuity,'-r', label='continuity')
	ax.semilogy(time,x_momentum,'-b', label='x-momentum')
	ax.semilogy(time,y_momentum,'-k', label='y-momentum')
	ax.semilogy(time,z_momentum,'-m', label='z-momentum')
	ax.semilogy(time,energy,'-g',     label='energy')
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles,labels)

def AppendToResidualsPlot(time, continuity, x_momentum, y_momentum, z_momentum, energy, ax):
	ax.semilogy(time,continuity,'-r')
	ax.semilogy(time,x_momentum,'-b')
	ax.semilogy(time,y_momentum,'-k')
	ax.semilogy(time,z_momentum,'-m')
	ax.semilogy(time,energy,'-g')


def UpdateResidualsPlot(fileName, skip_data, ax):
	N,time,continuity,x_momentum,y_momentum,z_momentum,energy = GatherNewResidualsValues(fileName, skip_data)

	if ( skip_data == 0 ):
		NewResidualsPlot(time, continuity, x_momentum, y_momentum, z_momentum, energy, ax)
	else:
		AppendToResidualsPlot(time, continuity, x_momentum, y_momentum, z_momentum, energy, ax)


	return (N)
