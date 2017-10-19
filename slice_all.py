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
base_dir = "/Users/morganmacleod/DATA/athenaruns/pm_envelope/smr_RL_hr_lr/" #"/home/morganmacleod/DATA/athenaruns/pm_envelope/convergence_studies/smr_dr_RL_0.2/"
m1 = 0.631686
m2 = 0.3
G=1

file_list = glob(base_dir+"HSE.out1.000[3-9][0-9].athdf")

mycm = plt.cm.bone_r

mylevel=2

vmin = -8

####################################


def get_midplane_theta(myfile,level=0):
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)

    # get closest to midplane value
    return dblank['x2v'][ np.argmin(np.abs(dblank['x2v']-np.pi/2.) ) ]



orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
t1=ou.get_t1(orb)
print "t1=",t1





x2slicevalue=get_midplane_theta(file_list[0],level=mylevel)
print "Slicing at x2=",x2slicevalue





for i,myfile in enumerate(file_list):
    print i, myfile 
    
    # read the data
    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                     x1_max=7.5,x2_min=x2slicevalue,x2_max=x2slicevalue)
    t = d['Time']
    
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
    
    theta_rot = -np.arctan2(y2,x2)
    #print y2,x2, theta_rot/(2.*np.pi)*360

    xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
    yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)

    plt.figure(figsize=(6,5))
    im=plt.pcolormesh(
               ou.get_plot_array_midplane(xrot),
               ou.get_plot_array_midplane(yrot),
               ou.get_plot_array_midplane(np.log10(d['rho'][:,len(d['x2v'])/2,:]) ),
              cmap=mycm,
              vmin=vmin,vmax=0,rasterized=True)
    
    plt.plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
         (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
         'w+')
    
    plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    plt.xlabel(r"$x'/R_1$")
    plt.ylabel(r"$y'/R_1$")
    plt.colorbar(im,label=r"$\log \ \rho \ \ [M_1/R_1^3]$")
    plt.savefig("snapshots/density_"+str(i)+".png",bbox_inches='tight',dpi=150)
    plt.clf()
