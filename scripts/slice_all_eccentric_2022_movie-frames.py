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
parser.add_argument("--filestart",help='int file number to start at',default=0,type=int)
parser.add_argument("--filestop",help='int file number to stop at',default=None,type=int)
parser.add_argument("--fileskip",help='int file number to skip',default=1,type=int)
parser.add_argument("--fn",help='str begining of filename',default='HSE.out1')
parser.add_argument("--level",help='int level to read file',default=None)
args = parser.parse_args()
base_dir=args.base_dir
filestart=args.filestart
filestop =args.filestop
fileskip =args.fileskip
fnstart = args.fn 
mylevel = args.level

G=6.674e-8
rsun = 6.955e10
file_list = sorted(glob(base_dir+fnstart+".[0-9][0-9][0-9][0-9][0-9].athdf"))
file_list = file_list[filestart:filestop:fileskip]
print (file_list)

#mylevel=None
####################################


def get_midplane_theta(myfile,level=0):
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)

    # get closest to midplane value
    return dblank['x2v'][ np.argmin(np.abs(dblank['x2v']-np.pi/2.) ) ]



orb = ou.read_trackfile(base_dir+"pm_trackfile.dat")


x2slicevalue=get_midplane_theta(file_list[0],level=mylevel)
print ("Slicing at x2=",x2slicevalue)



for i,myfile in enumerate(file_list):
    print (i, myfile )
    
    # read the data
    d = ou.read_data(myfile,orb,level=mylevel,get_cartesian=True,get_torque=False,get_energy=False,
                     x2_min=x2slicevalue,x2_max=x2slicevalue)
    t = d['Time']
    
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
    x2p_com = x2 - rcom[0]
    y2p_com = y2 - rcom[1]
    
    mycm = plt.cm.twilight
    fig,ax = plt.subplots(figsize=(9,6))
    im=plt.pcolormesh(
        ou.get_plot_array_midplane((d['x'][:,0,:]-rcom[0])/rsun),
        ou.get_plot_array_midplane((d['y'][:,0,:]-rcom[1])/rsun),
        ou.get_plot_array_midplane(np.log10(d['rho'][:,0,:]) ),
        cmap=mycm,
        vmin=-13,vmax=-3,rasterized=True)

    ax.add_patch(plt.Circle([(x2-rcom[0])/rsun,(y2-rcom[1])/rsun],radius=4,color='white') )
    

    #plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(-29,29),color='k',fontsize='small')
    plt.axis('equal')
    plt.xlabel(r"$x  \ \ [R_\odot]$")
    plt.ylabel(r"$y  \ \ [R_\odot]$")
    #plt.xlim(-lim,lim)
    #plt.ylim(-lim,lim)
    ax.set_ylim(-125,125)
    ax.set_xlim(-125,200)
    plt.colorbar(im,label=r"$\log_{10} \left( \rho \right)$   [g cm$^{-3}$]",extend='both')
    plt.savefig(myfile[:-6]+"_density_slice_com.png",bbox_inches='tight',dpi=100)
    plt.close()


    lim = 60
    fig,ax = plt.subplots(figsize=(9,7))
    im=plt.pcolormesh(
        ou.get_plot_array_midplane((d['x'][:,0,:])/rsun),
        ou.get_plot_array_midplane((d['y'][:,0,:])/rsun),
        ou.get_plot_array_midplane(np.log10(d['rho'][:,0,:]) ),
        cmap=mycm,
        vmin=-13,vmax=-3,rasterized=True)

    ax.add_patch(plt.Circle([(x2)/rsun,(y2)/rsun],radius=4,color='white') )
    

    #plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(-29,29),color='k',fontsize='small')
    plt.axis('equal')
    plt.xlabel(r"$x'  \ \ [R_\odot]$")
    plt.ylabel(r"$y'  \ \ [R_\odot]$")
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim)
    #ax.set_ylim(-9e12/rsun,7e12/rsun)
    #ax.set_xlim(-9e12/rsun,1.4e13/rsun)
    plt.colorbar(im,label=r"$\log_{10} \left( \rho \right)$   [g cm$^{-3}$]",extend='both')
    plt.savefig(myfile[:-6]+"_density_slice_zoom.png",bbox_inches='tight',dpi=100)
    plt.close()

    
    
