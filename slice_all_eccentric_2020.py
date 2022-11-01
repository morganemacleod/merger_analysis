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

args = parser.parse_args()
base_dir=args.base_dir
output_dir=args.output_dir

G=1
file_list = sorted(glob(base_dir+"HSE.out4.[0-9][0-9][0-9][0-9][0-9].athdf"))
file_list = file_list[1665::]
print file_list

mylevel=2
####################################


def get_midplane_theta(myfile,level=0):
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)

    # get closest to midplane value
    return dblank['x2v'][ np.argmin(np.abs(dblank['x2v']-np.pi/2.) ) ]



orb = ou.read_trackfile(base_dir+"pm_trackfile.dat")


x2slicevalue=get_midplane_theta(file_list[0],level=mylevel)
print "Slicing at x2=",x2slicevalue



for i,myfile in enumerate(file_list):
    print i, myfile 
    
    # read the data
    d = ou.read_data(myfile,orb,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,get_energy=False,
                     x2_min=x2slicevalue,x2_max=x2slicevalue,profile_file=base_dir+"hse_profile.dat",gamma=1.35)
    t = d['Time']
    
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
    x2p_com = x2 - rcom[0]
    y2p_com = y2 - rcom[1]
    
    vmin=-10
    mycm = plt.cm.magma
    plt.figure(figsize=(12,9))
    im=plt.pcolormesh(
        ou.get_plot_array_midplane(d['x'][:,len(d['x2v'])/2,:]-rcom[0]),
        ou.get_plot_array_midplane(d['y'][:,len(d['x2v'])/2,:]-rcom[1]),
        ou.get_plot_array_midplane(np.log10(d['rho'][:,len(d['x2v'])/2,:]) ),
        cmap=mycm,
        vmin=vmin,vmax=0,rasterized=True)

    plt.plot(x2p_com,y2p_com,'w*',markersize=3)
    

    plt.annotate(r"$t=$"+str(np.round(t,decimals=1)),(-29,29),color='k',fontsize='small')
    plt.axis('equal')
    plt.xlabel(r"$x'/R_1$")
    plt.ylabel(r"$y'/R_1$")
    plt.xlim(-3,12)
    plt.ylim(-7.5,7.5)
    plt.colorbar(im,label=r"$\log_{10} \left( \rho \right)$")
    plt.savefig(output_dir+"density_midplane_"+str(i)+".png",bbox_inches='tight',dpi=100)
    plt.close()

