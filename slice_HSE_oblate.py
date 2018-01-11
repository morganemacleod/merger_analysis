import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
import OrbitAnalysisUtils as ou

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

  
im=plt.pcolormesh(x,z,vel1,
                  cmap=mycm,
                  vmin=-0.1,vmax=0.1,rasterized=True)

plt.contour(x,z,rho,
            levels=[1.e-5],colors='k',linewidths=0.5)


plt.colorbar(im,label='Radial Velocity $[(G M_1 / R_1)^{1/2}]$')

plt.axis('equal')
plt.xlim(-2,2)
plt.ylim(-2,2)


plt.xlabel(r"$x/R_1$")
plt.ylabel(r"$z/R_1$")

plt.savefig(output_dir+"hse_3D_vel.pdf",dpi=300,bbox_inches='tight')
plt.clf()


##############
# DENSITY
##############
mycm = plt.cm.magma

# VERTICAL

x,z,rho = ou.get_plot_array_vertical("rho",0,
                       myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                       G=1,rsoft2=0.05,level=0,x1_max=3)

im=plt.pcolormesh(x,z,np.log10(rho),
                  cmap=mycm,
                  vmin=-8,vmax=0,rasterized=True)


plt.colorbar(im,label=r'$\log \ \rho \ \ [M_1 / R_1^3]$')


plt.xlabel(r"$x/R_1$")
plt.ylabel(r"$z/R_1$")

plt.axis('equal')
plt.xlim(-2,2)
plt.ylim(-2,2)

plt.savefig(output_dir+"hse_3D_density.pdf",dpi=300,bbox_inches='tight')
plt.clf()


