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
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/syncRL/fcorot_series/a206_res24_fc10/"

output_dir = "paper_figures/"

m1 = 0.410103
m2 = 0.3
G=1

full_file_list = glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf")

file_list = [full_file_list[-322],
             full_file_list[-62],
             full_file_list[-32],
             full_file_list[-23],
             full_file_list[-13],
             full_file_list[-1]]




mylevel=2

#vars = ['density','pressure','entropy', 'torque','cs','h']
#vars = ['h']
vars = ['density','pressure','entropy', 'torque','cs','h']

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

    mycm = plt.cm.magma
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
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,
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
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='w',fontsize='small')

        grid[i].axis('equal')
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$y'$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log_{10} \left(\rho \right)$")


    print "\n\n saving .... " , output_dir+"density_timeseries_midplane.pdf \n\n"
    plt.savefig(output_dir+"density_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.close()

        




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
        # read the data
        x,z,val = ou.get_plot_array_vertical("rho",theta_rot,
                       myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                             G=1,rsoft2=0.05,level=0,x1_max=7.5)

        im=grid[i].pcolormesh(x-np.linalg.norm(rcom),z,np.log10(val),
                  cmap=mycm,
                  vmin=vmin,vmax=0,rasterized=True)

        grid[i].plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='w',fontsize='small')

        grid[i].axis('equal')
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$z$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log_{10}\left( \rho \right)$")


    print "\n\n saving ... ", output_dir+"density_timeseries_vertical.pdf \n\n"
    plt.savefig(output_dir+"density_timeseries_vertical.pdf",bbox_inches='tight',dpi=300)
    plt.close()



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
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=True,
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
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='w',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$y'$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log_{10}\left(P \right)$")

    print "\n\n saving ... ", output_dir+"pressure_timeseries_midplane.pdf\n\n"
    plt.savefig(output_dir+"pressure_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.close()

        




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
        x,z,val = ou.get_plot_array_vertical("press",theta_rot,
                                             myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                             G=1,rsoft2=0.05,level=0,x1_max=7.5)

        im=grid[i].pcolormesh(x-np.linalg.norm(rcom),z,np.log10(val),
                  cmap=mycm,
                  vmin=vmin,vmax=0,rasterized=True)
        
        grid[i].plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='w',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$z$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log_{10} \left( P  \right)$")
        
    print "\n\n saving... ",output_dir+"pressure_timeseries_vertical.pdf \n\n"
    plt.savefig(output_dir+"pressure_timeseries_vertical.pdf",bbox_inches='tight',dpi=300)
    plt.close()



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
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=True,
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
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$y'$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\ln \left(P/\rho^\gamma \right)$")
    

    print "\n\n saving ... ",output_dir+"entropy_timeseries_midplane.pdf \n\n"
    plt.savefig(output_dir+"entropy_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.close()

        




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
        x,z,rho = ou.get_plot_array_vertical("rho",theta_rot,
                                             myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                             G=1,rsoft2=0.05,level=0,x1_max=7.5)
        x,z,pres = ou.get_plot_array_vertical("rho",theta_rot,
                                             myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                             G=1,rsoft2=0.05,level=0,x1_max=7.5)
        val = pres/rho**(5./3.)
        
        im=grid[i].pcolormesh(x-np.linalg.norm(rcom),z,np.log(val),
                  cmap=mycm,
                  vmin=vmin,vmax=vmax,rasterized=True)
        
        grid[i].plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$z$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\ln \left(P/\rho^\gamma \right)$")
            
    print "\n\n saving ... ",output_dir+"entropy_timeseries_vertical.pdf \n\n"
    plt.savefig(output_dir+"entropy_timeseries_vertical.pdf",bbox_inches='tight',dpi=300)
    plt.close()





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
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=True,
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
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$y'$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log_{10} \left( c_s \right)$")
    

    print "\n\n saving ... ",output_dir+"soundspeed_timeseries_midplane.pdf \n\n"
    plt.savefig(output_dir+"soundspeed_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.close()

        




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
        x,z,rho = ou.get_plot_array_vertical("rho",theta_rot,
                       myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                             G=1,rsoft2=0.05,level=0,x1_max=7.5)
        x,z,pres = ou.get_plot_array_vertical("press",theta_rot,
                       myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                             G=1,rsoft2=0.05,level=0,x1_max=7.5)

        val = np.sqrt(5./3. * pres/rho)

        im=grid[i].pcolormesh(x-np.linalg.norm(rcom),z,np.log10(val),
                  cmap=mycm,
                  vmin=vmin,vmax=0,rasterized=True)
        
        grid[i].plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$z$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\log_{10} \left( c_s \right)$")
            
    print "\n\n saving ... ",output_dir+"soundspeed_timeseries_vertical.pdf \n\n"
    plt.savefig(output_dir+"soundspeed_timeseries_vertical.pdf",bbox_inches='tight',dpi=300)
    plt.close()


###############
# TORQUE
###############
if 'torque' in vars:
    print "   TORQUE   "
    mycm = plt.cm.PiYG

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
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=True,
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

        grid[i].plot(0,0,
                     'kx',markersize=3)
        
        grid[i].plot((0-rcom[0])*np.cos(theta_rot)-(0-rcom[1])*np.sin(theta_rot),
                     (0-rcom[0])*np.sin(theta_rot)+(0-rcom[1])*np.cos(theta_rot),
                     'b*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-1.8,2.),color='k',fontsize='small')
        
        grid[i].set_xlim(-2,3)
        grid[i].set_ylim(-2.5,2.5)
        grid[i].set_xticks([-2,0,2])
        grid[i].set_yticks([-2,0,2])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$y'$")

        fig.colorbar(im,cax=grid.cbar_axes[i],label=r'$\hat z$ torque density')


    print "\n\n saving ... ",output_dir+"torque_timeseries_midplane.pdf\n\n"
    plt.savefig(output_dir+"torque_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.close()


####################################
### Specific Angular Momentum
####################################

if 'h' in vars:

    print "   SPECIFIC ANGULAR MOMENTUM     "
    
    vmax = 3.0
    vmin = 0.0
    mycm = plt.cm.YlGnBu_r

    
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
        d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x2_min=x2slicevalue,x2_max=x2slicevalue)
        t = d['Time']
        
        rcom,vcom = ou.rcom_vcom(orb,t)
        x2,y2,z2 = ou.pos_secondary(orb,t)
        
        theta_rot = -np.arctan2(y2,x2)
        #print y2,x2, theta_rot/(2.*np.pi)*360
        
        xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
        yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)

        d['h'] = (d['x']-rcom[0])*(d['vy']-vcom[1]) - (d['y']-rcom[1])*(d['vx']-vcom[0])


        im=grid[i].pcolormesh(
            ou.get_plot_array_midplane(xrot),
            ou.get_plot_array_midplane(yrot),
            ou.get_plot_array_midplane(d['h'][:,len(d['x2v'])/2,:] ),
            cmap=mycm,
            vmin=vmin,vmax=vmax,rasterized=True)

        grid[i].contour(ou.get_plot_array_midplane(xrot),
                    ou.get_plot_array_midplane(yrot),
                    ou.get_plot_array_midplane(np.log10(d['rho'][:,len(d['x2v'])/2,:])) ,
                        levels=[-5,-4,-3,-2,-1],colors='k',linewidths=0.5,linestyles='solid')
    
        grid[i].plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
                     (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
                     'w*',markersize=3)
        
        grid[i].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-4,4),color='k',fontsize='small')
        
        grid[i].set_xlim(-5,5)
        grid[i].set_ylim(-5,5)
        grid[i].set_xticks([-4,0,4])
        grid[i].set_yticks([-4,0,4])
        grid[i].set_xlabel(r"$x'$")
        grid[i].set_ylabel(r"$y'$")
        cb = grid.cbar_axes[i].colorbar(im)
        cb.set_label_text(r"$\hat z$ specific angular momentum")
    

    print "\n\n saving ... ",output_dir+"specmom_timeseries_midplane.pdf \n\n"
    plt.savefig(output_dir+"specmom_timeseries_midplane.pdf",bbox_inches='tight',dpi=300)
    plt.close()
