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
import matplotlib.colors as colors

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


#file_list = glob(base_dir+"HSE.out1.0023[0-8].athdf")
file_list = [base_dir+"HSE.out1.00030.athdf",
             base_dir+"HSE.out1.00043.athdf",
             base_dir+"HSE.out1.00053.athdf",
            base_dir+"HSE.out1.00067.athdf", 
            base_dir+"HSE.out1.00077.athdf",
            base_dir+"HSE.out1.00081.athdf",
            base_dir+"HSE.out1.00083.athdf",
            base_dir+"HSE.out1.00084.athdf",
            base_dir+"HSE.out1.00085.athdf"]

mycm = plt.cm.PiYG

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

fig = plt.figure(1,figsize=(8,11))
nrows = 3
ncols = 3
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(nrows,ncols),  # creates 2x2 grid of axes
                 axes_pad=0.1,  # pad between axes in inch.
                 cbar_mode='single',
                 cbar_location='right',
                 cbar_size="2%",
                 cbar_pad="2%",)



for i,myfile in enumerate(file_list):
    print i, myfile 
    
    # read the data
    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=True,
                     x1_max=7.5,x2_min=x2slicevalue,x2_max=x2slicevalue)

    t = d['Time']
    
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
    
    theta_rot = -np.arctan2(y2,x2)
    #print y2,x2, theta_rot/(2.*np.pi)*360

    xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
    yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)


    im=grid[i].pcolormesh(
               ou.get_plot_array_midplane(xrot),
               ou.get_plot_array_midplane(yrot),
               ou.get_plot_array_midplane(d['torque_dens_1_z'][:,len(d['x2v'])/2,:] +
                                          d['torque_dens_2_z'][:,len(d['x2v'])/2,:]) ,
               cmap=mycm,rasterized=True,
               norm=colors.SymLogNorm(linthresh=1.e-6,vmin=-1.0, vmax=1.0))
    
    grid[i].plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
         (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
         'w+')
    
    grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
    
    grid[i].set_xlim(-5,5)
    grid[i].set_ylim(-5,5)
    grid[i].set_xticks([-4,0,4])
    grid[i].set_yticks([-4,0,4])
    grid[i].set_xlabel(r"$x'/R_1$")
    grid[i].set_ylabel(r"$y'/R_1$")

    fig.colorbar(im,cax=grid.cbar_axes[i],label='torque density (z) [UNITS]')

    
plt.savefig("paper_figures/torque_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
plt.clf()

