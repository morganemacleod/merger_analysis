import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
import OrbitAnalysisUtils as ou
from mpl_toolkits.axes_grid1 import ImageGrid

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




###### SIMULATION PARAMS   #########
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/hse/res24_rot/"

output_dir = "paper_figures/"

m1 = 0.410103
m2 = 0.3
G=1

myfile = base_dir+"HSE.out1.00011.athdf"



mylevel=2

vmin = -8

####################################

def get_midplane_theta(myfile,level=0):
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)

    # get closest to midplane value
    return dblank['x2v'][ np.argmin(np.abs(dblank['x2v']-np.pi/2.) ) ]




orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")


fig = plt.figure(1,figsize=(10,9))
nrows = 2
ncols = 1
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(nrows,ncols),  # creates 2x2 grid of axes
                 axes_pad=0.3,  # pad between axes in inch.
                 cbar_mode='each',
                 cbar_location='right',
                 cbar_size="6%",
                 cbar_pad="5%",)



##############
# RADIAL VELOCITY
##############
mycm = plt.cm.RdBu_r

x,z,vel1 = ou.get_plot_array_vertical("vel1",0,
                       myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                       G=1,rsoft2=0.05,level=0,x1_max=3)

x,z,rho = ou.get_plot_array_vertical("rho",0,
                       myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                       G=1,rsoft2=0.05,level=0,x1_max=3)

  
im=grid[1].pcolormesh(x,z,vel1,
                  cmap=mycm,
                  vmin=-0.1,vmax=0.1,rasterized=True)

grid[1].contour(x,z,rho,
            levels=[1.e-5],colors='k',linewidths=0.5)


cb=grid.cbar_axes[1].colorbar(im)
cb.set_label_text(r'radial velocity')

grid[1].axis('equal')
grid[1].set_xlim(-2,2)
grid[1].set_ylim(-2,2)


grid[1].set_xlabel(r"$x$")
grid[1].set_ylabel(r"$z$")

#plt.savefig(output_dir+"hse_3D_vel.pdf",dpi=300,bbox_inches='tight')
#plt.clf()


##############
# DENSITY
##############
mycm = plt.cm.magma

# VERTICAL

x,z,rho = ou.get_plot_array_vertical("rho",0,
                       myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                       G=1,rsoft2=0.05,level=0,x1_max=3)

im=grid[0].pcolormesh(x,z,np.log10(rho),
                  cmap=mycm,
                  vmin=-8,vmax=0,rasterized=True)


cb=grid.cbar_axes[0].colorbar(im)
cb.set_label_text(r'$\log_{10} \left( \rho \right)$')


grid[0].set_xlabel(r"$x$")
grid[0].set_ylabel(r"$z$")

grid[0].axis('equal')
grid[0].set_xlim(-2,2)
grid[0].set_ylim(-2,2)

plt.savefig(output_dir+"hse_3D.pdf",dpi=300,bbox_inches='tight')
plt.clf()


