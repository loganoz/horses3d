# HORSES3D Python Real-Time Monitors

This is a PYTHON utility to plot real time monitors. This package requires *matplotlib* in the environment `PATH`. The most popular `python` open-source distribution with *matplotlib* is Anaconda:

https://www.anaconda.com/download

## Usage

Monitors are invoked by calling the function `plotMonitors`, in the `Monitors.py` module. This can be done either using a small python script or directly from the python command line. This function needs a set of monitor file names from a list. The whole procedure is:

1. Load the Monitors environment. This can be done by sourcing the `configure-python.sh` script:
```bash
$ source HORSES_TOP_LEVEL/utils/PythonUtilities/configure-python
```
2. `cd` to the current case folder
3. open a python shell by entering
```bash
$ python
```
4. And, insider the python shell:

  1. Import the monitors module:
```PYTHON
   >> import Monitors
```
  2. Set the files to plot:
```PYTHON
   >> fileNames = ["file.residuals","file.monitor.volume"]
```
  3. Call the `plotMonitors` function
```PYTHON
   Monitors.plotMonitors(fileNames)
```

  A python `example.py` script will look like:
```PYTHON
   import Monitors
   fileNames = ["file.residuals","file.monitor.volume"]
   Monitors.plotMonitors(fileNames)
```

  and can be executed using:
```bash
$ python example.py
```
