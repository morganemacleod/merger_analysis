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
base_dir =  "/Volumes/DATAVolume/athenaruns/pm_envelope/smr_RL_hr_lr2_f15/"
m1 = 0.631686
m2 = 0.3
G=1

file_list = glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf")

output_dir = "snapshots/sma2_f15/"
mycm = plt.cm.bone_r

mylevel=2

vmin = -8

vars = ['rho','press','cs','etot','torque','entropy']
#vars = ['entropy']

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
    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=True,get_energy=True,
                     x1_max=8.,x2_min=x2slicevalue,x2_max=x2slicevalue,profile_file=base_dir+"hse_profile.dat")
    t = d['Time']
    
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
    
    theta_rot = -np.arctan2(y2,x2)
    #print y2,x2, theta_rot/(2.*np.pi)*360

    xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
    yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)

    x2p_com = (x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot)
    y2p_com = (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot)


    if 'rho' in vars:
        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(np.log10(d['rho'][:,len(d['x2v'])/2,:]) ),
            cmap=mycm,
            vmin=vmin,vmax=0,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"$\log \ \rho \ \ [M_1/R_1^3]$")
        plt.savefig(output_dir+"density_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()

        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(np.log10(d['rho'][:,len(d['x2v'])/2,:]) ),
            cmap=mycm,
            vmin=vmin,vmax=0,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(x2p_com-1,x2p_com+1)
        plt.ylim(y2p_com-1,y2p_com+1)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"$\log \ \rho \ \ [M_1/R_1^3]$")
        plt.savefig(output_dir+"density_zoom_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()


    
    if 'press' in vars:
        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(np.log10(d['press'][:,len(d['x2v'])/2,:]) ),
            cmap=plt.cm.viridis,
            vmin=vmin,vmax=0,rasterized=True)

        plt.plot(x2p_com,y2p_com,'w*')

        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"$\log \ P$")
        plt.savefig(output_dir+"pressure_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()

        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(np.log10(d['press'][:,len(d['x2v'])/2,:]) ),
            cmap=plt.cm.viridis,
            vmin=vmin,vmax=0,rasterized=True)

        plt.plot(x2p_com,y2p_com,'w*')

        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(x2p_com-1,x2p_com+1)
        plt.ylim(y2p_com-1,y2p_com+1)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"$\log \ P$")
        plt.savefig(output_dir+"pressure_zoom_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()


    if 'cs' in vars:
        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(np.log10( np.sqrt(5./3.* d['press'][:,len(d['x2v'])/2,:] /
                                          d['rho'][:,len(d['x2v'])/2,:])  ) ),
            cmap=plt.cm.plasma,
            vmin=-2,vmax=0.3,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"$\log c_s [(GM_1/R_1)^{1/2}]$")
        plt.savefig(output_dir+"sound_speed_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()

        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(np.log10( np.sqrt(5./3.* d['press'][:,len(d['x2v'])/2,:] /
                                          d['rho'][:,len(d['x2v'])/2,:])  ) ),
            cmap=plt.cm.plasma,
            vmin=-2,vmax=0.3,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(x2p_com-1,x2p_com+1)
        plt.ylim(y2p_com-1,y2p_com+1)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"$\log c_s [(GM_1/R_1)^{1/2}]$")
        plt.savefig(output_dir+"sound_speed_zoom_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()


    if 'entropy' in vars:
        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(np.log( d['press'][:,len(d['x2v'])/2,:] / d['rho'][:,len(d['x2v'])/2,:]**(5./3.)   ) ),
            cmap=plt.cm.RdYlBu_r,
            vmin=-0.25,vmax=5.5,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"entropy $\ln \left(P/\rho^\gamma \right)$")
        plt.savefig(output_dir+"entropy_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()

        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(np.log( d['press'][:,len(d['x2v'])/2,:] / d['rho'][:,len(d['x2v'])/2,:]**(5./3.)   ) ),
            cmap=plt.cm.RdYlBu_r,
            vmin=-0.25,vmax=5.5,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(x2p_com-1,x2p_com+1)
        plt.ylim(y2p_com-1,y2p_com+1)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"entropy $\ln \left(P/\rho^\gamma \right)$")
        plt.savefig(output_dir+"entropy_zoom_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()

    
    if 'etot' in vars:
        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(d['etot'][:,len(d['x2v'])/2,:] /
                                          d['rho'][:,len(d['x2v'])/2,:]) ,
            cmap=plt.cm.RdBu_r,rasterized=True,
            vmin=-1,vmax=1)

        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"specific binding energy $\log \ e_{\rm tot} \ \ [GM_1/R_1]$")
        plt.savefig(output_dir+"specificBE_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()

        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(d['etot'][:,len(d['x2v'])/2,:] /
                                          d['rho'][:,len(d['x2v'])/2,:]) ,
            cmap=plt.cm.RdBu_r,rasterized=True,
            vmin=-1,vmax=1)

        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(x2p_com-1,x2p_com+1)
        plt.ylim(y2p_com-1,y2p_com+1)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"specific binding energy $\log \ e_{\rm tot} \ \ [GM_1/R_1]$")
        plt.savefig(output_dir+"specificBE_zoom_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()


    if 'torque' in vars:
        mynorm =colors.SymLogNorm(linthresh=1.e-6,vmin=-1.0, vmax=1.0)
        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
               ou.get_plot_array_midplane(xrot),
               ou.get_plot_array_midplane(yrot),
               ou.get_plot_array_midplane(d['torque_dens_1_z'][:,len(d['x2v'])/2,:] +
                                          d['torque_dens_2_z'][:,len(d['x2v'])/2,:]) ,
               cmap=plt.cm.PiYG,rasterized=True,
               norm=mynorm)

        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"torque density (z)")
        plt.savefig(output_dir+"torque_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()


        plt.figure(figsize=(6.2,5))
        im=plt.pcolormesh(
               ou.get_plot_array_midplane(xrot),
               ou.get_plot_array_midplane(yrot),
               ou.get_plot_array_midplane(d['torque_dens_1_z'][:,len(d['x2v'])/2,:] +
                                          d['torque_dens_2_z'][:,len(d['x2v'])/2,:]) ,
               cmap=plt.cm.PiYG,rasterized=True,
               norm=mynorm)

        
        plt.plot(x2p_com,y2p_com,'w*')
    
        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        plt.xlim(x2p_com-1,x2p_com+1)
        plt.ylim(y2p_com-1,y2p_com+1)
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.colorbar(im,label=r"torque density (z)")
        plt.savefig(output_dir+"torque_zoom_"+str(i)+".png",bbox_inches='tight',dpi=300)
        plt.close()
