import numpy as np
import matplotlib.pyplot as plt
from merger_analysis import athena_read as ar
from merger_analysis import OrbitAnalysisUtils as ou
from Constants import Constants
#import seaborn as sns
#import deepdish as dd
#from astropy.table import Table
#from glob import glob
#from mpl_toolkits.axes_grid1 import ImageGrid
#from tqdm.auto import tqdm

c=Constants()

#%matplotlib inline

plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.2
plt.rcParams['legend.labelspacing'] = 0.2
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 16



