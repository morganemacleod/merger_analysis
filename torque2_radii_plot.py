import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
import OrbitAnalysisUtils as ou
import argparse

# set some global options
plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 16



### SIMULATION PARAMS ######
m1 = 0.410103
m2 = 0.3

base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/syncRL/fcorot_series/a206_res24_fc10/"

############################

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")



###
# MASS IN RADII
###
tr = ascii.read("torque2_radii_time.dat")
tr['sep'] = np.interp(tr['time'],orb['time'],orb['sep'])


#plt.figure(figsize=(8,5))

mycm = LinearSegmentedColormap.from_list("mycm", ["SandyBrown","Black"])
labels=[r"$r_2<0.1$",r"$r_2<0.2$",r"$r_2<0.3$",r"$r_2<0.4$",r"$r_2<0.6$"]

cols = np.array(tr.colnames)
cols = cols[[1,2,3,4,6]]
print "plotting columns,",cols

for i,colname in enumerate(cols):
    plt.plot(tr['sep'],
             tr[colname],
            color=mycm(float(i)/len(cols) ),
            label=labels[i],lw=2)

plt.xlabel("separation")
plt.grid()


# TOTAL GRAV TORQUE ON THE TWO COMPONENTS
orb['t2'] = np.cross(orb['r']-orb['rcom'],m2*orb['agas2'])[:,2]

plt.plot(orb['sep'],orb['t2'],'k',label=r'total',lw=2)

plt.ylabel(r'$\tau_{\rm grav,2}$')
plt.legend(loc=0,frameon=True)

plt.savefig("paper_figures/torque_2_spatial.pdf",bbox_inches='tight')
