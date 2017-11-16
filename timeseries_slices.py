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
base_dir = "/Users/morganmacleod/DATA/athenaruns/pm_envelope/smr_RL_hr_lr2/" 

output_dir = "paper_figures/"

m1 = 0.631686
m2 = 0.3
G=1

full_file_list = glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf")

file_list = [full_file_list[-100],
             full_file_list[-35],
             full_file_list[-15],
             full_file_list[-5],
             full_file_list[-2],
             full_file_list[-1]]




mylevel=2

vars = ['density','pressure','entropy', 'torque','cs']
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


####################################
### DENSITY
####################################

if 'density' in vars:

    print " MAKING DENSITY TIMESERIES "

    mycm = plt.cm.bone_r
    vmin = -8

    
    fig = plt.figure(1,figsize=(10,9))
    nrows = 2
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
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
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
            ou.get_plot_array_midplane(np.log10(d['rho'][:,len(d['x2v'])/2,:]) ),
            cmap=mycm,
            vmin=vmin,vmax=0,rasterized=True)
    
        grid[i].plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
                     (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$y'/R_1$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log \ \rho \ \ [M_1/R_1^3]$")


    print "\n\n saving .... " , output_dir+"density_timeseries_midplane.pdf \n\n"
    plt.savefig(output_dir+"density_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.clf()

        




    fig = plt.figure(1,figsize=(10,9))
    nrows = 2
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
        
        dblank=ar.athdf(myfile,level=mylevel,quantities=[],subsample=True)
        t=dblank['Time']
        rcom,vcom = ou.rcom_vcom(orb,t)
        x2,y2,z2 = ou.pos_secondary(orb,t)
        
        theta_rot = -np.arctan2(y2,x2)
        x3slicevalue = dblank['x3v'][np.argmin(np.abs(dblank['x3v']+theta_rot))]
        print "x3=",x3slicevalue
        # read the data
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        
        im=grid[i].pcolormesh(d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                              d['z'][0,:,:],
                              np.log10(d['rho'][0,:,:]),
                              cmap=mycm,
                              vmin=vmin,vmax=0,rasterized=True)
        
        if(x3slicevalue<0):
            x3slicevalue += np.pi
        else:
            x3slicevalue -= np.pi
        print "x3=",x3slicevalue
        # read the data
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        im=grid[i].pcolormesh(-d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                            d['z'][0,:,:],
                              np.log10(d['rho'][0,:,:]),
                              cmap=mycm,
                              vmin=vmin,vmax=0,rasterized=True)
        
        grid[i].plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$z/R_1$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log \ \rho \ \ [M_1/R_1^3]$")


    print "\n\n saving ... ", output_dir+"density_timeseries_vertical.pdf \n\n"
    plt.savefig(output_dir+"density_timeseries_vertical.pdf",bbox_inches='tight',dpi=300)
    plt.clf()



####################################
### PRESSURE
####################################

if 'pressure' in vars:

    print " PRESSURE PLOT  "

    mycm = plt.cm.viridis
    vmin = -8
    
    fig = plt.figure(1,figsize=(10,9))
    nrows = 2
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
            ou.get_plot_array_midplane(np.log10(d['press'][:,len(d['x2v'])/2,:]) ),
            cmap=mycm,
            vmin=vmin,vmax=0,rasterized=True)
    
        grid[i].plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
                     (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$y'/R_1$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log \ \rho \ \ [M_1/R_1^3]$")

    print "\n\n saving ... ", output_dir+"pressure_timeseries_midplane.pdf\n\n"
    plt.savefig(output_dir+"pressure_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.clf()

        




    fig = plt.figure(1,figsize=(10,9))
    nrows = 2
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
            
        dblank=ar.athdf(myfile,level=mylevel,quantities=[],subsample=True)
        t=dblank['Time']
        rcom,vcom = ou.rcom_vcom(orb,t)
        x2,y2,z2 = ou.pos_secondary(orb,t)
        
        theta_rot = -np.arctan2(y2,x2)
        x3slicevalue = dblank['x3v'][np.argmin(np.abs(dblank['x3v']+theta_rot))]
        print "x3=",x3slicevalue
        # read the data
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        
        im=grid[i].pcolormesh(d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                              d['z'][0,:,:],
                              np.log10(d['press'][0,:,:]),
                              cmap=mycm,
                              vmin=vmin,vmax=0,rasterized=True)
        
        if(x3slicevalue<0):
            x3slicevalue += np.pi
        else:
            x3slicevalue -= np.pi
        print "x3=",x3slicevalue
        # read the data
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        im=grid[i].pcolormesh(-d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                              d['z'][0,:,:],
                              np.log10(d['press'][0,:,:]),
                              cmap=mycm,
                              vmin=vmin,vmax=0,rasterized=True)
        
        grid[i].plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$z/R_1$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log \ \rho \ \ [M_1/R_1^3]$")
        
    print "\n\n saving... ",output_dir+"pressure_timeseries_vertical.pdf \n\n"
    plt.savefig(output_dir+"pressure_timeseries_vertical.pdf",bbox_inches='tight',dpi=300)
    plt.clf()



####################################
### ENTROPY
####################################

if 'entropy' in vars:

    print "   ENTROPY     "
    
    vmax = 10.0
    vmin = -0.25
    mycm = plt.cm.RdYlBu_r

    
    fig = plt.figure(1,figsize=(10,9))
    nrows = 2
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
            ou.get_plot_array_midplane(np.log( d['press'][:,len(d['x2v'])/2,:] / d['rho'][:,len(d['x2v'])/2,:]**(5./3.)   ) ),
            cmap=mycm,
            vmin=vmin,vmax=vmax,rasterized=True)
    
        grid[i].plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
                     (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$y'/R_1$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log \ \rho \ \ [M_1/R_1^3]$")
    

    print "\n\n saving ... ",output_dir+"entropy_timeseries_midplane.pdf \n\n"
    plt.savefig(output_dir+"entropy_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.clf()

        




    fig = plt.figure(1,figsize=(10,9))
    nrows = 2
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
        
        dblank=ar.athdf(myfile,level=mylevel,quantities=[],subsample=True)
        t=dblank['Time']
        rcom,vcom = ou.rcom_vcom(orb,t)
        x2,y2,z2 = ou.pos_secondary(orb,t)
        
        theta_rot = -np.arctan2(y2,x2)
        x3slicevalue = dblank['x3v'][np.argmin(np.abs(dblank['x3v']+theta_rot))]
        print "x3=",x3slicevalue
        # read the data
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        
        im=grid[i].pcolormesh(d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                              d['z'][0,:,:],
                              np.log(d['press'][0,:,:]/d['rho'][0,:,:]**(5./3.) ),
                              cmap=mycm,
                              vmin=vmin,vmax=vmax,rasterized=True)
        
        if(x3slicevalue<0):
            x3slicevalue += np.pi
        else:
            x3slicevalue -= np.pi
        print "x3=",x3slicevalue
        # read the data
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        im=grid[i].pcolormesh(-d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                              d['z'][0,:,:],
                              np.log(d['press'][0,:,:]/d['rho'][0,:,:]**(5./3.) ),
                              cmap=mycm,
                              vmin=vmin,vmax=vmax,rasterized=True)
        
        grid[i].plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$z/R_1$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log \ \rho \ \ [M_1/R_1^3]$")
            
    print "\n\n saving ... ",output_dir+"entropy_timeseries_vertical.pdf \n\n"
    plt.savefig(output_dir+"entropy_timeseries_vertical.pdf",bbox_inches='tight',dpi=300)
    plt.clf()





####################################
### SOUND SPEED
####################################

if 'cs' in vars:

    print "   SOUND SPEED     "
    
    vmax = 0.3
    vmin = -2
    mycm = plt.cm.plasma

    
    fig = plt.figure(1,figsize=(10,9))
    nrows = 2
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
            ou.get_plot_array_midplane(np.log10((5./3.)* np.sqrt( d['press'][:,len(d['x2v'])/2,:] / d['rho'][:,len(d['x2v'])/2,:]  )) ),
            cmap=mycm,
            vmin=vmin,vmax=vmax,rasterized=True)
    
        grid[i].plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
                     (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$y'/R_1$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log \ c_s \ \ [(GM_1/R_1)^{1/2}]$")
    

    print "\n\n saving ... ",output_dir+"soundspeed_timeseries_midplane.pdf \n\n"
    plt.savefig(output_dir+"soundspeed_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.clf()

        




    fig = plt.figure(1,figsize=(10,9))
    nrows = 2
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
        
        dblank=ar.athdf(myfile,level=mylevel,quantities=[],subsample=True)
        t=dblank['Time']
        rcom,vcom = ou.rcom_vcom(orb,t)
        x2,y2,z2 = ou.pos_secondary(orb,t)
        
        theta_rot = -np.arctan2(y2,x2)
        x3slicevalue = dblank['x3v'][np.argmin(np.abs(dblank['x3v']+theta_rot))]
        print "x3=",x3slicevalue
        # read the data
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        
        im=grid[i].pcolormesh(d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                              d['z'][0,:,:],
                              np.log10(np.sqrt(5./3.*d['press'][0,:,:]/d['rho'][0,:,:] )),
                              cmap=mycm,
                              vmin=vmin,vmax=0,rasterized=True)
        
        if(x3slicevalue<0):
            x3slicevalue += np.pi
        else:
            x3slicevalue -= np.pi
        print "x3=",x3slicevalue
        # read the data
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        im=grid[i].pcolormesh(-d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                              d['z'][0,:,:],
                              np.log10(np.sqrt(5./3.*d['press'][0,:,:]/d['rho'][0,:,:] )),
                              cmap=mycm,
                              vmin=vmin,vmax=vmax,rasterized=True)
        
        grid[i].plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$z/R_1$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log \ c_s \ \ [(GM_1/R_1)^{1/2}]$")
            
    print "\n\n saving ... ",output_dir+"soundspeed_timeseries_vertical.pdf \n\n"
    plt.savefig(output_dir+"soundspeed_timeseries_vertical.pdf",bbox_inches='tight',dpi=300)
    plt.clf()



if 'torque' in vars:
    print "   TORQUE   "
    mycm = plt.cm.PiYG
    vmin = -8
    fig = plt.figure(1,figsize=(8,11))
    nrows = 2
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
                     'w*',markersize=3)
    
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,3.5),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'/R_1$")
        grid[i].set_ylabel(r"$y'/R_1$")

        fig.colorbar(im,cax=grid.cbar_axes[i],label='torque density (z) [UNITS]')


    print "\n\n saving ... ",output_dir+"torque_timeseries_midplane.pdf\n\n"
    plt.savefig(output_dir+"torque_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.clf()
