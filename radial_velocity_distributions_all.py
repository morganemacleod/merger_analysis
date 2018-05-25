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
import seaborn as sns

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
parser = argparse.ArgumentParser(description='Read m1,m2, input/output directories')

parser.add_argument("m1",type=float,help="mass of particle m1")
parser.add_argument("m2",type=float,help="mass of particle m2")

parser.add_argument("--base_dir", help="data directory (should end with / )")
parser.add_argument("--output_dir", help="directory to save figures/output (should end with / )")

args = parser.parse_args()
m1=args.m1
m2=args.m2
base_dir=args.base_dir
output_dir=args.output_dir

G=1
file_list = sorted(glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf"))
file_list = file_list[30::]
print file_list

mylevel=0

####################################


orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
t1=ou.get_t1(orb)
print "t1=",t1


for i,myfile in enumerate(file_list):
    print i, myfile 

    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=0,get_cartesian=True,get_torque=False,
                 get_energy=True,    
                 x1_max=30)
    t = d['Time']

    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
    
    theta_rot = -np.arctan2(y2,x2)
    
    sma = np.interp(t,orb['time'],orb['sep'])
    Omega = np.interp(t,orb['time'],orb['vmag'])/sma 
    print sma,Omega

    # ROTATE POSITIONS, VELOCITIES TO COROTATING FRAME
    xrot = (d['x']-rcom[0])*np.cos(theta_rot) - (d['y']-rcom[1])*np.sin(theta_rot)
    yrot = (d['x']-rcom[0])*np.sin(theta_rot) + (d['y']-rcom[1])*np.cos(theta_rot)

    vxrot = (d['vx']-vcom[0])*np.cos(theta_rot) - (d['vy']-vcom[1])*np.sin(theta_rot)
    vyrot = (d['vx']-vcom[0])*np.sin(theta_rot) + (d['vy']-vcom[1])*np.cos(theta_rot)
    

    phicom = np.arctan2(yrot,xrot)

    vphi = - Omega*np.sqrt(xrot**2 + yrot**2)
    vxrotC = vxrot - vphi*np.sin(phicom)
    vyrotC = vyrot + vphi*np.cos(phicom)

    x2rot = (x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot)
    y2rot = (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot)

    EJ = (0.5*(vxrotC**2 +vyrotC**2 + d['vz']**2)
          +((d['epotp']+d['epotg'])/d['rho'])
          -0.5*((-Omega*yrot)**2 + (Omega*xrot)**2)
    )

    hz = xrot*vyrot - yrot*vxrot
    
    d['vr_com'] = np.sqrt(vxrot**2 + vyrot**2 +d['vz']**2)



    # VELOCITY PLOTS
    mycm= sns.cubehelix_palette(start=0.9, rot=-1,light=0.95,dark=.1,as_cmap=True)
    rl = np.linspace(1,30,100)
    lc = 'C4'


    ######################
    #   Velocity
    ######################

    histx =  np.sqrt(xrot**2 + yrot**2 +d['z']**2).flatten()  
    histy = d['vr_com'].flatten()
    dm = (d['dvol']*d['rho']).flatten()
    
    plt.hist2d(histx,histy,weights=dm,
               bins=60,norm=colors.LogNorm(vmin=1.e-8,vmax=1.e-2),
               range=[[1,30],[0,1]],
    cmap=mycm)
    plt.colorbar(label="mass per zone",ticks=10**np.linspace(-12,-2,11))
    
    plt.plot(rl,np.sqrt(2*1.3/rl),lw=2,color=lc)
    plt.annotate(r"$v_{esc}(r)$",(25,0.35),color=lc)
    
    plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(0.7,0.92),xycoords='axes fraction',
                 fontsize='small',color='k',backgroundcolor='white')

    plt.xlabel(r"$ r_{\rm com}/R_1 $")
    plt.ylabel(r"$v_{r,com}$")
    plt.savefig(output_dir+"velocity_dist_full_"+str(i)+".png",bbox_inches='tight',dpi=150)
    plt.close()



    ## THREE PANEL LINEAR ##
    fig=plt.figure(figsize=(5.5,12))

    nbins=60

    midplane = (np.abs(d['gx2v']-np.pi/2)< 2*np.pi* 30/360.).flatten()
    intermediate = ((np.abs(d['gx2v']-np.pi/2) >  2*np.pi* 30/360.) &
                    (np.abs(d['gx2v']-np.pi/2) <= 2*np.pi* 60/360.)).flatten()
    pole = (np.abs(d['gx2v']-np.pi/2) > 2*np.pi* 60/360.).flatten()

    
    histx =  np.sqrt(xrot**2 + yrot**2 +d['z']**2).flatten()  
    histy = d['vr_com'].flatten()
    dm = (d['dvol']*d['rho']).flatten()



    ### MIDPLANE
    plt.subplot(311)
    plt.hist2d(histx[midplane],histy[midplane],weights=dm[midplane],
               bins=nbins,norm=colors.LogNorm(vmin=1.e-8,vmax=1.e-2),
               range=[[1,30],[0,1]],
               cmap=mycm)
    

    
    plt.plot(rl,np.sqrt(2*1.3/rl),lw=2,color=lc)
    plt.annotate(r"$v_{esc}(r)$",(25,0.35),color=lc)
    
    plt.annotate(r"lat$<30^\degree$",(0.35,0.92),xycoords='axes fraction')
    
    plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(0.7,0.92),xycoords='axes fraction',
                 fontsize='small',color='k')

    #plt.xlabel(r"$ r_{\rm com}/R_1 $")
    plt.ylabel(r"$v_{r,com}$")
    plt.xticks(visible=False)
    plt.yticks([0.2,0.4,0.6,0.8,1.0])

    ### INTERMEDIATE
    plt.subplot(312)
    plt.hist2d(histx[intermediate],histy[intermediate],weights=dm[intermediate],
               bins=nbins,norm=colors.LogNorm(vmin=1.e-8,vmax=1.e-2),
               range=[[1,30],[0,1]],
               cmap=mycm)
    
    plt.plot(rl,np.sqrt(2*1.3/rl),lw=2,color=lc)
    plt.annotate(r"$v_{esc}(r)$",(25,0.35),color=lc)
    
    plt.annotate(r"$30^\degree<$lat$<60^\degree$",(0.35,0.92),xycoords='axes fraction')
    
    #plt.xlabel(r"$ r_{\rm com}/R_1 $")
    plt.ylabel(r"$v_{r,com}$")
    plt.xticks(visible=False)
    plt.yticks([0.2,0.4,0.6,0.8,1.0])
    

    ### POLE
    plt.subplot(313)
    c,xe,ye,im=plt.hist2d(histx[pole],histy[pole],weights=dm[pole],
                          bins=nbins,norm=colors.LogNorm(vmin=1.e-8,vmax=1.e-2),
                          range=[[1,30],[0,1]],
                          cmap=mycm)



    plt.plot(rl,np.sqrt(2*1.3/rl),lw=2,color=lc)
    plt.annotate(r"$v_{esc}(r)$",(25,0.35),color=lc)
    
    plt.annotate(r"lat$>60^\degree$",(0.35,0.92),xycoords='axes fraction')
    
    plt.xlabel(r"$ r_{\rm com}/R_1 $")
    plt.ylabel(r"$v_{r,com}$")
    

    fig.subplots_adjust(hspace=0.0, right=0.9)
    cbar_ax = fig.add_axes([0.925, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label="mass per zone")

    plt.savefig(output_dir+"velocity_dist_split_linear_"+str(i)+".png",bbox_inches='tight',dpi=150)
    plt.close()





    ### LOGARITHMIC ###

    fig=plt.figure(figsize=(5.5,12))

    nbins=60
    myrange = [[0,1.5],[-1,0.5]]
    
    midplane = (np.abs(d['gx2v']-np.pi/2)< 2*np.pi* 30/360.).flatten()
    intermediate = ((np.abs(d['gx2v']-np.pi/2) >  2*np.pi* 30/360.) &
                    (np.abs(d['gx2v']-np.pi/2) <= 2*np.pi* 60/360.)).flatten()
    pole = (np.abs(d['gx2v']-np.pi/2) > 2*np.pi* 60/360.).flatten()
    
    histx =  np.log10(np.sqrt(xrot**2 + yrot**2 +d['z']**2)).flatten()  
    histy = np.log10(d['vr_com']).flatten()
    dm = (d['dvol']*d['rho']).flatten()



    ### MIDPLANE
    plt.subplot(311)
    plt.hist2d(histx[midplane],histy[midplane],weights=dm[midplane],
               bins=nbins,norm=colors.LogNorm(vmin=1.e-8,vmax=1.e-2),
               range=myrange,
               cmap=mycm)
    
    
    
    plt.plot(np.log10(rl),np.log10(np.sqrt(2*1.3/rl) ),lw=2,color=lc)
    plt.annotate(r"$v_{esc}(r)$",(1.1,-0.28),color=lc)
    
    plt.annotate(r"lat$<30^\degree$",(0.35,0.92),xycoords='axes fraction')
    
    plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(0.7,0.92),xycoords='axes fraction',
                 fontsize='small',color='k')
    
    #plt.xlabel(r"$ r_{\rm com}/R_1 $")
    plt.ylabel(r"$\log_{10 }(v_{r,com})$")
    plt.xticks(visible=False)
    plt.yticks([-0.75,-0.5,-0.25,0,0.25,0.5])
    
    ### INTERMEDIATE
    plt.subplot(312)
    plt.hist2d(histx[intermediate],histy[intermediate],weights=dm[intermediate],
               bins=nbins,norm=colors.LogNorm(vmin=1.e-8,vmax=1.e-2),
               range=myrange,
               cmap=mycm)

    plt.plot(np.log10(rl),np.log10(np.sqrt(2*1.3/rl) ),lw=2,color=lc)
    plt.annotate(r"$v_{esc}(r)$",(1.1,-0.28),color=lc)
    
    plt.annotate(r"$30^\degree<$lat$<60^\degree$",(0.35,0.92),xycoords='axes fraction')
    
    #plt.xlabel(r"$ r_{\rm com}/R_1 $")
    plt.ylabel(r"$\log_{10 }(v_{r,com})$")
    plt.xticks(visible=False)
    plt.yticks([-0.75,-0.5,-0.25,0,0.25,0.5])
    
    
    ### POLE
    plt.subplot(313)
    c,xe,ye,im=plt.hist2d(histx[pole],histy[pole],weights=dm[pole],
                          bins=nbins,norm=colors.LogNorm(vmin=1.e-8,vmax=1.e-2),
                          range=myrange,
                          cmap=mycm)
    
    

    plt.plot(np.log10(rl),np.log10(np.sqrt(2*1.3/rl) ),lw=2,color=lc)
    plt.annotate(r"$v_{esc}(r)$",(1.1,-0.28),color=lc)
    
    plt.annotate(r"lat$>60^\degree$",(0.35,0.92),xycoords='axes fraction')
    
    plt.xlabel(r"$\log_{10 } ( r_{\rm com}/R_1 )$")
    plt.ylabel(r"$v_{r,com}$")
    
    
    fig.subplots_adjust(hspace=0.0, right=0.9)
    cbar_ax = fig.add_axes([0.925, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label="mass per zone")
    plt.savefig(output_dir+"velocity_dist_split_log_"+str(i)+".png",bbox_inches='tight',dpi=150)
    plt.close()

    



