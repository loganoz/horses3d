"""
This script generates 2D plots from Paraview-generated slice ASCII files in *.csv format.
Usage: <python> plot_paraview_slices.py infile_path
"""

import os
import json
import string
import argparse
import numpy as np
from typing import Union
from scipy.interpolate import griddata

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
   # "font.sans-serif": "Helvetica"
    "legend.handlelength": 1.0
})

from matplotlib import colors

def plot_field_value(data_path: list[str], nx: int, ny: int, 
                     xlabel: str, ylabel: Union[str, list[str]], xticks: list[float], yticks: list[float], xtick_labels: list[str], ytick_labels: list[str],
                     axes: Union[list[int], list[list[int]]], field: Union[int,list[int]], field_label: str, figwidth: float, out_path: str, 
                     plotarray : list[int] = [0,0], y_to_x_ratio : float = 1.0, invert_xaxis = False,
                     cmap : str = 'coolwarm', cmin : Union[None, float] = None, cmax : Union[None, float] = None, cticks : Union[None, list[float]] = None, cnorm : str = 'linear',
                     fontsize : float = 9, figlabel_loc : Union[None, list[float]] = None, ymask : Union[None, list[str]] = None, xmask : Union[None, list[str]] = None,
                     xlim : Union[None, list[float], list[list[float]]] = None, ylim : Union[None, list[float], list[list[float]]] = None) -> None:
    
    N = len(data_path)
    datasets = []
    xmin_list = []
    xmax_list = []
    ymin_list = []
    ymax_list = []

    plt.rcParams.update({"font.size": fontsize})

    # Load datasets
    for i, path in enumerate(data_path):
        rawdata = np.loadtxt(path, skiprows=1, delimiter=',')
        if isinstance(axes[0], list):
            xdata = axes[0][i]
            ydata = axes[1][i]
        else:
            xdata = axes[0]
            ydata = axes[1]
        if isinstance(field, list):
            fielddata = field[i]
        else:
            fielddata = field
        x = rawdata[:, xdata]
        y = rawdata[:, ydata]
        field_value = rawdata[:, fielddata]
       
        # Get min and max values for representation
        if xlim is None:
            xmax = max(x)
            xmin = min(x)
        elif isinstance(xlim[0], list):
            xmax = min(xlim[i][1], max(x))
            xmin = max(xlim[i][0], min(x))
        else:   
            xmax = min(xlim[1], max(x))
            xmin = max(xlim[0], min(x))

        if ylim is None:
            ymax = max(y)
            ymin = min(y)
        elif isinstance(ylim[0], list):
            ymax = min(ylim[i][1], max(y))
            ymin = max(ylim[i][0], min(y))
        else:   
            ymax = min(ylim[1], max(y))
            ymin = max(ylim[0], min(y))

        # Define grid for interpolation
        grid_x, grid_y = np.mgrid[xmin:xmax:complex(nx), ymin:ymax:complex(ny)]

        # ################################ TMP ###################################
        if cnorm == 'log':
            field_value[field_value < 0] = 1e-8
        # ########################################################################

        # Interpolate the field values on the grid
        datasets.append(griddata((x, y), field_value, (grid_x, grid_y), method='linear'))

        # ################################ TMP ###################################
        # datasets[-1][datasets[-1] < 0] = 1e-16
        # ########################################################################

        # Apply mask (if defined)
        if ymask is not None:
            curve = eval(ymask[i].replace('x', 'grid_x'))
            mask_points = grid_y < curve
            datasets[i][mask_points] = np.nan
        if xmask is not None:
            curve = eval(xmask[i].replace('y', 'grid_y'))
            mask_points = grid_x < curve
            datasets[i][mask_points] = np.nan

        # Store min and max values for x and y
        xmin_list.append(xmin)
        xmax_list.append(xmax)
        ymin_list.append(ymin)
        ymax_list.append(ymax)

    sharey = False
    if plotarray != [0,0]:
        nrows = plotarray[0]
        ncols = plotarray[1]
        if ncols > 1:
            sharey = True
    else:
        nrows = N
        ncols = 1

    r_row = [1.0] * nrows
    irow = 0
    for i in range(nrows*ncols):
        if i % ncols == 0:
            r_row[irow] = (ymax_list[i]-ymin_list[i])/(xmax_list[i]-xmin_list[i]) * y_to_x_ratio
            irow += 1
    #r = (max(y)-min(y))/(max(x)-min(x)) * y_to_x_ratio
    fraction = 0.025 # fraction of the sum of x axes width that will be used for the colorbar
    textheight = fontsize/72 # in inches, 1pt is 1/72 inches -> 10pt
    figheight = sum(r_row) * (figwidth - 4*textheight) / (ncols*(1+fraction)) + 2*textheight
    #figheight = r * nrows * (figwidth - 4*textheight) / (ncols*(1+fraction)) + 2*textheight
    aspect = sum(r_row) / fraction / ncols

    if cmin is not None:
        vmin = cmin
    else:
        vmin = np.nanmin(datasets)

    if cmax is not None:
        vmax = cmax
    else:
        vmax = np.nanmax(datasets)

    if cnorm == 'linear':
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
    if cnorm == 'log':
        norm=colors.LogNorm(vmin=vmin, vmax=vmax)
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=sharey, figsize=(figwidth, figheight), gridspec_kw={'height_ratios': r_row}, layout='compressed')
    images = []
    if N == 1:
        axs = np.array([axs])
    
    count = 0
    for ax, data in zip(axs.flat, datasets):
        images.append(ax.imshow(data.T, extent=(xmin_list[count], xmax_list[count], ymin_list[count], ymax_list[count]), origin='lower', cmap=cmap, aspect='auto', norm=norm))
        count += 1
    
    if ncols == 1 and nrows >= 1:
        axs[-1].set_xticks(ticks=xticks, labels=xtick_labels)
        axs[-1].set_xlabel(xlabel)
        count = 0
        for i in range(N):
            if isinstance(yticks[0], list) and isinstance(ytick_labels[0], list):
                axs[i].set_yticks(ticks=yticks[count], labels=ytick_labels[count])
            else:
                axs[i].set_yticks(ticks=yticks, labels=ytick_labels)
            if isinstance(ylabel, list):
                axs[i].set_ylabel(ylabel[count])
            else:
                axs[i].set_ylabel(ylabel)
            count += 1
    elif nrows == 1 and ncols > 1:
        axs[0].set_yticks(ticks=yticks, labels=ytick_labels)
        axs[0].set_ylabel(ylabel)
        for i in range(N):
            axs[i].set_xticks(ticks=xticks, labels=xtick_labels)
            axs[i].set_xlabel(xlabel)
    else:
        for icol in range(ncols):
            axs[-1, icol].set_xticks(ticks=xticks, labels=xtick_labels)
            axs[-1, icol].set_xlabel(xlabel)
        for irow in range(nrows):
            axs[irow, 0].set_yticks(ticks=yticks, labels=ytick_labels)
            axs[irow, 0].set_ylabel(ylabel)
    
    if figlabel_loc:
        labels = string.ascii_lowercase   
        for i, ax in enumerate(axs.flat):
            ax.text(figlabel_loc[0], figlabel_loc[1], '('+labels[i]+')', transform=ax.transAxes, fontsize=fontsize, verticalalignment='top')
    
    fig.align_ylabels(axs)
    fig.align_xlabels(axs)    
    # fig.supylabel(ylabel)

    fig.colorbar(images[-1], ax=axs.ravel().tolist(), fraction=fraction, aspect=aspect, label=field_label, pad=0.015, ticks=cticks)
    
    if invert_xaxis:
        for ax in axs.flat:
            ax.invert_xaxis()
    
    plt.savefig(out_path, dpi=300)

def main():

    parser = argparse.ArgumentParser(description="Generates 2D plots from Paraview-generated slice ASCII files in *.csv format.")
    
    # Define postional arguments
    parser.add_argument('infile_path', type=str, help='Path to the input *.json file.')

    # Parse arguments
    args = parser.parse_args()

    # Change the working directory to that where the input JSON file is located
    os.chdir(os.path.dirname(os.path.abspath(args.infile_path)))

    # Load json file
    with open(os.path.abspath(args.infile_path), 'r') as f:
        inputs_dict = json.load(f)

    for key in inputs_dict:

        fig_dict = inputs_dict[key]
        
        kwargs = {}
        if 'plotArray' in fig_dict:
            kwargs['plotarray'] = fig_dict['plotArray']
        if 'yxRatio' in fig_dict:
            kwargs['y_to_x_ratio'] = fig_dict['yxRatio']
        if 'invertXaxis' in fig_dict:
            kwargs['invert_xaxis'] = fig_dict['invertXaxis']
        if 'colorMap' in fig_dict:
            kwargs['cmap'] = fig_dict['colorMap']
        if 'colorMapMin' in fig_dict:
            kwargs['cmin'] = fig_dict['colorMapMin']
        if 'colorMapMax' in fig_dict:
            kwargs['cmax'] = fig_dict['colorMapMax']
        if 'colorBarTicks' in fig_dict:
            kwargs['cticks'] = fig_dict['colorBarTicks']
        if 'colorMapNorm' in fig_dict:
            kwargs['cnorm'] = fig_dict['colorMapNorm']
        if 'fontsize' in fig_dict:
            kwargs['fontsize'] = fig_dict['fontsize']
        if 'figLabelLoc' in fig_dict:
            kwargs['figlabel_loc'] = fig_dict['figLabelLoc']
        if 'xLim' in fig_dict:
            kwargs['xlim'] = fig_dict['xLim']
        if 'yLim' in fig_dict:
            kwargs['ylim'] = fig_dict['yLim']
        if 'xMask' in fig_dict:
            kwargs['xmask'] = fig_dict['xMask']
        if 'yMask' in fig_dict:
            kwargs['ymask'] = fig_dict['yMask']
        

        plot_field_value(fig_dict['dataPathList'], fig_dict['nx'], fig_dict['ny'], fig_dict['xlabel'], fig_dict['ylabel'], fig_dict['xticks'], fig_dict['yticks'], fig_dict['xtickLabels'], fig_dict['ytickLabels'], 
                         [fig_dict['xCol'], fig_dict['yCol']], fig_dict['dataCol'], fig_dict['dataLabel'], fig_dict['figwidth'], fig_dict['outputPath'], **kwargs)

if __name__ == "__main__":
    main()