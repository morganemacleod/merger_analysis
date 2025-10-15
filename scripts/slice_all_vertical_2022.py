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
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.colors as colors
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



###### SIMULATION PARAMS   #########
parser = argparse.ArgumentParser(description='Read input/output directories')

parser.add_argument("--base_dir", help="data directory (should end with / )")
parser.add_argument("--output_dir", help="directory to save figures/output (should end with / )")
parser.add_argument("--rstar_init",help="initial stellar polar radius",type=float)
parser.add_argument("--filestart",help='int file number to start at',default=0,type=int)
parser.add_argument("--filestop",help='int file number to stop at',default=-1,type=int)
parser.add_argument("--fileskip",help='int file number to skip',default=1,type=int)


args = parser.parse_args()
base_dir=args.base_dir
output_dir=args.output_dir
rstar_init = args.rstar_init
filestart=args.filestart
filestop =args.filestop
fileskip =args.fileskip


G=6.674e-8
file_list = sorted(glob(base_dir+"HSE.out1.[0-9][0-9][0-9][0-9][0-9].athdf"))
#file_list = file_list[466::]
file_list = file_list[filestart:filestop:fileskip]
print (file_list)

mylevel=None
####################################


def get_midplane_theta(myfile,level=0):
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)

    # get closest to midplane value
    return dblank['x2v'][ np.argmin(np.abs(dblank['x2v']-np.pi/2.) ) ]

def get_phi_pi(myfile,level=0):
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)

    # get closest to midplane value
    return dblank['x3v'][ np.argmin(np.abs(dblank['x3v']-np.pi) ) ]



orb = ou.read_trackfile(base_dir+"pm_trackfile.dat")


x3slicevalue=get_phi_pi(file_list[0],level=mylevel)
print ("Slicing at x3=",x3slicevalue)



for i,myfile in enumerate(file_list):
    print (i, myfile )
    
    # read the data
    d = ou.read_data(myfile,orb,level=mylevel,get_cartesian=True,get_torque=False,get_energy=False,
                     x3_min=x3slicevalue,x3_max=x3slicevalue,x1_max=6.e12)
    t = d['Time']
    
    plt.figure(figsize=(6,6))
    plt.pcolormesh(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                   d['z'][0,:,:],
                   np.log10((d['rho'])[0,:,:]),cmap='magma',rasterized=True)
    plt.colorbar()
    plt.contour(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                d['z'][0,:,:],
                np.log10((d['rho'])[0,:,:]),levels=[-10,-8,-6,-4,-3],colors='k',linestyles='-')
    plt.axis('equal')
    plt.xlim(0,3e12)
    plt.ylim(-2.e12,2e12)
    plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(1e11,1.8e12),color='k',fontsize='small')
    plt.savefig(output_dir+"vertical_density_"+"%05d"%i+".png",bbox_inches='tight',dpi=100)
    plt.close()


    plt.figure(figsize=(6,6))
    plt.pcolormesh(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                   d['z'][0,:,:],
                   np.log10(np.abs(d['vel1'])/np.sqrt(1.35*d['press']/d['rho']))[0,:,:],vmin=-3,vmax=3,cmap='RdBu' ,rasterized=True)
    plt.colorbar()
    plt.contour(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                d['z'][0,:,:],
                np.log10((d['rho'])[0,:,:]),levels=[-10,-8,-6,-4,-3],colors='k',linestyles='-')
    plt.axis('equal')
    plt.xlim(0,3e12)
    plt.ylim(-2.e12,2e12)
    plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(1e11,1.8e12),color='k',fontsize='small')
    plt.savefig(output_dir+"vertical_mach_"+"%05d"%i+".png",bbox_inches='tight',dpi=100)
    plt.close()

    plt.figure(figsize=(6,6))
    plt.pcolormesh(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                   d['z'][0,:,:],
                   d['vel3'][0,:,:],cmap='Spectral' )
    plt.colorbar()
    plt.contour(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
            d['z'][0,:,:],
            np.log10((d['rho'])[0,:,:]),levels=[-10,-8,-6,-4,-3],colors='k',linestyles='-')
    plt.axis('equal')
    plt.xlim(0,3e12)
    plt.ylim(-2.e12,2e12)
    plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(1e11,1.8e12),color='k',fontsize='small')
    plt.savefig(output_dir+"vertical_vphi_"+"%05d"%i+".png",bbox_inches='tight',dpi=100)
    plt.close()


    plt.figure(figsize=(6,6))
    plt.pcolormesh(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                   d['z'][0,:,:],
                   (d['vel3']/np.sqrt(d['x']**2 + d['y']**2))[0,:,:]/1.45299e-05,cmap='Spectral',vmin=0,vmax=2 ,rasterized=True)
    plt.colorbar()
    plt.contour(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
            d['z'][0,:,:],
            np.log10((d['rho'])[0,:,:]),levels=[-10,-8,-6,-4,-3],colors='k',linestyles='-')
    plt.axis('equal')
    plt.xlim(0,3e12)
    plt.ylim(-2.e12,2e12)
    plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(1e11,1.8e12),color='k',fontsize='small')
    plt.savefig(output_dir+"vertical_omeganorm_"+"%05d"%i+".png",bbox_inches='tight',dpi=100)
    plt.close()

    plt.figure(figsize=(6,6))
    plt.pcolormesh(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                   d['z'][0,:,:],
                   (d['vel3']/np.sqrt(2.0/3.0*6.67e-8*6.905e+34/rstar_init))[0,:,:],cmap='Spectral',vmin=0,vmax=1 ,rasterized=True)
    plt.colorbar()
    plt.contour(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
            d['z'][0,:,:],
            np.log10((d['rho'])[0,:,:]),levels=[-10,-8,-6,-4,-3],colors='k',linestyles='-')
    plt.axis('equal')
    plt.xlim(0,3e12)
    plt.ylim(-2.e12,2e12)
    plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(1e11,1.8e12),color='k',fontsize='small')
    plt.savefig(output_dir+"vertical_vcrit_"+"%05d"%i+".png",bbox_inches='tight',dpi=100)
    plt.close()

    plt.figure(figsize=(6,6))
    plt.pcolormesh(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                   d['z'][0,:,:],
                   d['r0'][0,:,:],cmap='Spectral',vmin=0,vmax=1 ,rasterized=True)
    plt.colorbar()
    plt.contour(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
            d['z'][0,:,:],
            np.log10((d['rho'])[0,:,:]),levels=[-10,-8,-6,-4,-3],colors='k',linestyles='-')
    plt.axis('equal')
    plt.xlim(0,3e12)
    plt.ylim(-2.e12,2e12)
    plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(1e11,1.8e12),color='k',fontsize='small')
    plt.savefig(output_dir+"vertical_r0_"+"%05d"%i+".png",bbox_inches='tight',dpi=100)
    plt.close()

    plt.figure(figsize=(6,6))
    plt.pcolormesh(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
                   d['z'][0,:,:],
                   np.log10(np.sqrt(d['press']/d['rho'])/1.e5)[0,:,:],cmap='Spectral',vmin=1,vmax=3,rasterized=True)
    plt.colorbar()
    plt.contour(np.sqrt(d['x']**2 + d['y']**2)[0,:,:],
            d['z'][0,:,:],
            np.log10((d['rho'])[0,:,:]),levels=[-10,-8,-6,-4,-3],colors='k',linestyles='-')
    plt.axis('equal')
    plt.xlim(0,3e12)
    plt.ylim(-2.e12,2e12)
    plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(1e11,1.8e12),color='k',fontsize='small')
    plt.savefig(output_dir+"vertical_cs_"+"%05d"%i+".png",bbox_inches='tight',dpi=100)
    plt.close()
