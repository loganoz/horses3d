#!/usr/bin/env python
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from BenchUtilities import *

(elements,poly,L2error,elapsedTime,CPUTime) = GatherResults()
PlotBench(elements,poly,L2error,elapsedTime,CPUTime)

f = plt.gcf()
f.set_size_inches(11.69,8.27) 	#A4 paper size
plt.savefig("./RESULTS/BenchResults.pdf")


#mng = plt.get_current_fig_manager()
#mng.resize(*mng.window.maxsize())
plt.show()

