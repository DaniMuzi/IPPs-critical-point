#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt


import matplotlib.ticker as tkr
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

import matplotlib.cm as cm


from matplotlib import colors
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as patches

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)

# plt.rc('font',**{'family':'serif','serif':['Times']})          # Everyday life
plt.rc('font',**{'family':'serif','serif':['Helvetica']})    # Saving figures
plt.rc('text', usetex=True)

plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['axes.titlesize'] = 30
plt.rcParams['axes.labelsize'] = 35

plt.rcParams['legend.loc'] = 'best'
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.shadow'] = False
plt.rcParams['legend.framealpha'] = 0.0
plt.rcParams['legend.markerscale'] = 1.5
plt.rcParams['legend.handletextpad'] = 0.5

plt.rcParams['xtick.major.size'] = 9
plt.rcParams['ytick.major.size'] = 9
plt.rcParams['xtick.minor.size'] = 0
plt.rcParams['ytick.minor.size'] = 0

plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

import warnings
warnings.filterwarnings('ignore')
