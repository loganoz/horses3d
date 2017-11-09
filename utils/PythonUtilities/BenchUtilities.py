#!/usr/bin/env python
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def GatherResults():
	fid = open("./RESULTS/BenchResults.out",'r')
	elements = []
	poly = []
	L2error = []
	elapsedTime = []
	CPUTime = []
	for line in fid:
		raw_data = line.split()
		elements.append(int(raw_data[0]))
		poly.append(int(raw_data[1]))
		L2error.append(float(raw_data[2]))
		elapsedTime.append(float(raw_data[3]))
		CPUTime.append(float(raw_data[4]))

	return (elements,poly,L2error,elapsedTime,CPUTime)

def AppendL2Plot(ax,NDOF,poly,L2,element):
#
#	Add the plot
#	------------
	name = "$" + str(round((element)**(1/3))) + "^3$"
	L2plot = ax.plot(NDOF,L2,'-s',markersize=12,markerfacecolor=[1,1,1], mew=1.5,label=name)

	for i in range(0,len(NDOF)):
		ax.annotate(str(poly[i]),(NDOF[i],L2[i]),ha="center",va="center",size=8,weight="bold")

	return L2plot
	
def AppendTimePlot(axE,axC,NDOF,poly,elapsed,CPU,maxElapsed,element):
#
#	Add the plot
#	------------
	elapsedPlot = axE.plot(NDOF,elapsed/maxElapsed,'-s',markersize=12,markerfacecolor=[1,1,1], mew=1.5)
	for i in range(0,len(NDOF)):
		axE.annotate(str(poly[i]),(NDOF[i],elapsed[i]/maxElapsed),ha="center",va="center",size=8,weight="bold")
#
#	Plot the parallel efficiency
#	----------------------------
	quotient = np.array(CPU,dtype=np.float) / np.array(elapsed,dtype=np.float)

	type(elapsedPlot)
	type(elapsedPlot[0])
	color = elapsedPlot[0].get_color()
	axC.plot(NDOF,quotient,'--o',color=color,markersize=5,markerfacecolor=[1,1,1], mew=1.5)
	axC.set_xscale('log')
	axC.set_ylim([1,10])
	axC.set_xticks([])




def PlotL2error(ax,elements,poly,L2error):
#
#	Create plot
#	-----------
	ax.clear()
	currentElement = elements[0]
	L2 = []
	NDOF = []
	plot_poly = []
	for i in range(0,len(elements)):
		if ( currentElement != elements[i] ):
			AppendL2Plot(ax,NDOF,plot_poly,L2,currentElement)
			currentElement = elements[i]
			NDOF = []
			L2 = []
			plot_poly = []
		NDOF.append((poly[i]+1.0)**3.0*currentElement)
		L2.append(L2error[i])
		plot_poly.append(poly[i])
	AppendL2Plot(ax,NDOF,plot_poly,L2,currentElement)
	ax.set_xlabel('NDOF')
	ax.set_ylabel('$L^2$ error')
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.legend(loc='upper right',ncol=4,fancybox=True, shadow=True,bbox_to_anchor=(1.4,1.1))
	

def PlotTime(axE,axC,elements,poly,elapsedTime,CPUTime):
#
#	Compute maximum elapsed time
#	----------------------------
	maxTime = np.amax(elapsedTime)	
#
#	Create plot
#	-----------
	plt.figure(1)
	plt.rcParams['axes.linewidth'] = 1.5
	axE.clear()
	axC.clear()
	currentElement = elements[0]
	eT = []
	CPUT = []
	NDOF = []
	plot_poly = []
	for i in range(0,len(elements)):
		if ( currentElement != elements[i] ):
			AppendTimePlot(axE,axC,NDOF,plot_poly,eT,CPUT,maxTime,currentElement)
			currentElement = elements[i]
			NDOF = []
			eT = []
			CPUT = []
			plot_poly = []
		NDOF.append((poly[i]+1.0)**3.0*currentElement)
		eT.append(elapsedTime[i])
		CPUT.append(CPUTime[i])
		plot_poly.append(poly[i])
	AppendTimePlot(axE,axC,NDOF,plot_poly,eT,CPUT,maxTime,currentElement)
	label = "Maximum time: " + str(round(maxTime,2))
	axE.annotate(label,(1000,1))
	axE.set_xlabel('NDOF')
	axE.set_ylabel('Elapsed time / max(Elapsed time)')
	axE.yaxis.set_label_position("right")
	axE.set_yscale('log')
	axE.set_xscale('log')
	axC.set_ylabel('$\eta_{PARALLEL}$ (%)')
	axC.yaxis.set_label_position("right")

def PlotBench(elements,poly,L2,elapsedTime,CPUTime):
	fig = plt.figure(1)
#
#	Define the layout
#	-----------------
	gs = gridspec.GridSpec(3,2)
#
#	Plot the L2 error
#	-----------------
	ax1 = plt.subplot(gs[:, 0])
	PlotL2error(ax1,elements,poly,L2)
#
#	Plot the elapsed and CPU time
#	-----------------------------
	ax2 = plt.subplot(gs[1:3,1])
	ax3 = plt.subplot(gs[0,1])
	PlotTime(ax2,ax3,elements,poly,elapsedTime,CPUTime)
