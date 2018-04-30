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

mylevel=2


do_vertical = True
plot_vars = ['rho','entropy','mach','vr','divV']
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
    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,get_energy=False,
                     x2_min=x2slicevalue,x2_max=x2slicevalue,profile_file=base_dir+"hse_profile.dat")
    t = d['Time']
    
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
    x2p_com = x2 - rcom[0]
    y2p_com = y2 - rcom[1]


    if 'rho' in plot_vars:
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
    

        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
        plt.axis('equal')
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.xlim(-20,20)
        plt.ylim(-20,20)
        plt.colorbar(im,label=r"$\log_{10} \left( \rho \right)$")
        plt.savefig(output_dir+"density_midplane_"+str(i)+".png",bbox_inches='tight',dpi=100)
        plt.close()

    if 'entropy' in plot_vars:
        mycm = plt.cm.RdYlBu_r
        plt.figure(figsize=(12,9))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(d['x'][:,len(d['x2v'])/2,:]-rcom[0]),
            ou.get_plot_array_midplane(d['y'][:,len(d['x2v'])/2,:]-rcom[1]),
            ou.get_plot_array_midplane(np.log(d['press'][:,0,:]/d['rho'][:,0,:]**(5./3.) ) ),
            cmap=mycm,
            vmin=-1,vmax=10,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*',markersize=3)
    

        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
        plt.axis('equal')
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.xlim(-20,20)
        plt.ylim(-20,20)
        plt.colorbar(im,label=r"$\ln \left( P / \rho^\gamma \right)$")
        plt.savefig(output_dir+"entropy_midplane_"+str(i)+".png",bbox_inches='tight',dpi=100)
        plt.close()

    if 'mach' in plot_vars:
        mycm = plt.cm.RdBu
        cs = np.sqrt((5./3.)*d['press'][:,0,:]/d['rho'][:,0,:] )
        vmag = np.sqrt((d['vx'][:,0,:]-vcom[0])**2 + (d['vy'][:,0,:]-vcom[1])**2 + (d['vz'][:,0,:]-vcom[2])**2 )
        plt.figure(figsize=(12,9))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(d['x'][:,len(d['x2v'])/2,:]-rcom[0]),
            ou.get_plot_array_midplane(d['y'][:,len(d['x2v'])/2,:]-rcom[1]),
            ou.get_plot_array_midplane(np.log10(vmag/cs ) ),
            cmap=mycm,
            vmin=-2,vmax=2,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*',markersize=3)
    

        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
        plt.axis('equal')
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.xlim(-20,20)
        plt.ylim(-20,20)
        plt.colorbar(im,label=r"log mach number (inertial)")
        plt.savefig(output_dir+"mach_midplane_"+str(i)+".png",bbox_inches='tight',dpi=100)
        plt.close()

    if 'vr' in plot_vars:
        mycm = plt.cm.RdBu
        radius_com = np.sqrt((d['x']-rcom[0])**2 + (d['y']-rcom[1])**2 + (d['z']-rcom[2])**2)
        d['vr_com'] = ((d['vx']-vcom[0])*(d['x']-rcom[0])/radius_com +
                       (d['vy']-vcom[1])*(d['y']-rcom[1])/radius_com +
                       (d['vz']-vcom[2])*(d['z']-rcom[2])/radius_com )
  
        plt.figure(figsize=(12,9))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(d['x'][:,len(d['x2v'])/2,:]-rcom[0]),
            ou.get_plot_array_midplane(d['y'][:,len(d['x2v'])/2,:]-rcom[1]),
            ou.get_plot_array_midplane(d['vr_com'][:,0,:] ) ,
            cmap=mycm,
            vmin=-1,vmax=1,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*',markersize=3)
    

        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
        plt.axis('equal')
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.xlim(-20,20)
        plt.ylim(-20,20)
        plt.colorbar(im,label=r"radial velocity (inertial)")
        plt.savefig(output_dir+"vr_midplane_"+str(i)+".png",bbox_inches='tight',dpi=100)
        plt.close()

    if 'divV' in plot_vars:
        mycm = plt.cm.PRGn

        # r, phi terms here... 
        divV = ( d['gx1v']**(-2) *  np.gradient(d['gx1v']**2 * d['vel1'],axis=2)/np.gradient(d['gx1v'],axis=2) +
                 #1/(d['gx1v']*np.sin(d['gx2v'])) *  np.gradient(np.sin(d['gx2v']) *d['vel2'],axis=1)/np.gradient(d['gx2v'],axis=1) +
                 1/(d['gx1v']*np.sin(d['gx2v'])) *  np.gradient(d['vel3'],axis=0)/np.gradient(d['gx3v'],axis=0) )

        plt.figure(figsize=(12,9))
        im=plt.pcolormesh(
            ou.get_plot_array_midplane(d['x'][:,len(d['x2v'])/2,:]-rcom[0]),
            ou.get_plot_array_midplane(d['y'][:,len(d['x2v'])/2,:]-rcom[1]),
            ou.get_plot_array_midplane(divV[:,0,:] ) ,
            cmap=mycm,
            vmin=-0.5,vmax=0.5,rasterized=True)
        
        plt.plot(x2p_com,y2p_com,'w*',markersize=3)
    

        plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
        plt.axis('equal')
        plt.xlabel(r"$x'/R_1$")
        plt.ylabel(r"$y'/R_1$")
        plt.xlim(-20,20)
        plt.ylim(-20,20)
        plt.colorbar(im,label=r"$\nabla v$")
        plt.savefig(output_dir+"divV_midplane_"+str(i)+".png",bbox_inches='tight',dpi=100)
        plt.close()

    

    if do_vertical:
        # VERTICAL
        dblank=ar.athdf(myfile,level=mylevel,quantities=[],subsample=True)
        t=dblank['Time']
        rcom,vcom = ou.rcom_vcom(orb,t)
        x2,y2,z2 = ou.pos_secondary(orb,t)
        
        theta_rot = -np.arctan2(y2,x2)
        # read the data
        d={}
        x,z,d['rho'] = ou.get_plot_array_vertical("rho",theta_rot,
                                                  myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                  G=1,rsoft2=0.05,level=mylevel)
        x,z,d['x'] = ou.get_plot_array_vertical("x",theta_rot,
                                                myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                G=1,rsoft2=0.05,level=mylevel)
        x,z,d['y'] = ou.get_plot_array_vertical("y",theta_rot,
                                                myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                G=1,rsoft2=0.05,level=mylevel)
        x,z,d['z'] = ou.get_plot_array_vertical("z",theta_rot,
                                                myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                G=1,rsoft2=0.05,level=mylevel)
        
        x,z,d['press'] = ou.get_plot_array_vertical("press",theta_rot,
                                                    myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                    G=1,rsoft2=0.05,level=mylevel)
        
        x,z,d['vx'] = ou.get_plot_array_vertical("vx",theta_rot,
                                                 myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                 G=1,rsoft2=0.05,level=mylevel)
        
        x,z,d['vy'] = ou.get_plot_array_vertical("vy",theta_rot,
                                                 myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                 G=1,rsoft2=0.05,level=mylevel)
        
        
        x,z,d['vz'] = ou.get_plot_array_vertical("vz",theta_rot,
                                                 myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                 G=1,rsoft2=0.05,level=mylevel)

        x,z,d['gx1v'] = ou.get_plot_array_vertical("gx1v",theta_rot,
                                                   myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                   G=1,rsoft2=0.05,level=mylevel)
        x,z,d['gx2v'] = ou.get_plot_array_vertical("gx2v",theta_rot,
                                                   myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                   G=1,rsoft2=0.05,level=mylevel)
        
        x,z,d['vel1'] = ou.get_plot_array_vertical("vel1",theta_rot,
                                                   myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                   G=1,rsoft2=0.05,level=mylevel)
        x,z,d['vel2'] = ou.get_plot_array_vertical("vel2",theta_rot,
                                                   myfile,base_dir+"hse_profile.dat",orb,m1,m2,
                                                   G=1,rsoft2=0.05,level=mylevel)
        

        r = np.sqrt(x**2 + z**2)
        radius_com = np.sqrt((d['x']-rcom[0])**2 + (d['y']-rcom[1])**2 + (d['z']-rcom[2])**2)
        d['vr_com'] = ((d['vx']-vcom[0])*(d['x']-rcom[0])/radius_com +
                       (d['vy']-vcom[1])*(d['y']-rcom[1])/radius_com +
                       (d['vz']-vcom[2])*(d['z']-rcom[2])/radius_com )
        d['vmag_com'] = np.sqrt((d['vx']-vcom[0])**2 + (d['vy']-vcom[1])**2 + (d['vz']-vcom[2])**2 )
        



        if 'rho' in plot_vars:
            plt.figure(figsize=(12,9) )
            im=plt.pcolormesh(x-np.linalg.norm(rcom),z,np.log10(d['rho']),
                              cmap=plt.cm.magma,
                              rasterized=True,
                              vmin=-10,vmax=0)
        
            plt.plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
            
            plt.axis('equal')
            plt.xlabel(r"$x'/R_1$")
            plt.ylabel(r"$z/R_1$")
            plt.colorbar(im,label=r"$\log_{10} \left( \rho \right)$")
            
            plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
            
            plt.savefig(output_dir+"density_vert_"+str(i)+".png",bbox_inches='tight',dpi=100)
            plt.close()


        
        if 'entropy' in plot_vars:
            plt.figure(figsize=(12,9) )
            im=plt.pcolormesh(x-np.linalg.norm(rcom),z,np.log(d['press'] /d['rho']**(5./3.) ),
                              cmap=plt.cm.RdYlBu_r,
                              rasterized=True,
                              vmin=-1,vmax=10)
            
            plt.plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
            
            plt.axis('equal')
            plt.xlabel(r"$x'/R_1$")
            plt.ylabel(r"$z/R_1$")
            plt.colorbar(im,label=r"$\ln \left( P/ \rho^\gamma \right)$")
            
            plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
            
            plt.savefig(output_dir+"entropy_vert_"+str(i)+".png",bbox_inches='tight',dpi=100)
            plt.close()

        if 'mach' in plot_vars:
            plt.figure(figsize=(12,9) )
            im=plt.pcolormesh(x-np.linalg.norm(rcom),z,
                              np.log10( d['vmag_com']/np.sqrt((5./3.)*d['press'] /d['rho'] ) ),
                              cmap=plt.cm.RdBu,
                              rasterized=True,
                              vmin=-2,vmax=2)
        
            plt.plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
            
            plt.axis('equal')
            plt.xlabel(r"$x'/R_1$")
            plt.ylabel(r"$z/R_1$")
            plt.colorbar(im,label=r"log mach number (inertial)")
            
            plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
            
            plt.savefig(output_dir+"mach_vert_"+str(i)+".png",bbox_inches='tight',dpi=100)
            plt.close()


        if 'vr' in plot_vars:
            plt.figure(figsize=(12,9) )
            im=plt.pcolormesh(x-np.linalg.norm(rcom),z,
                              d['vr_com'],
                              cmap=plt.cm.RdBu,
                              rasterized=True,
                              vmin=-1,vmax=1)
        
            plt.plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
            
            plt.axis('equal')
            plt.xlabel(r"$x'/R_1$")
            plt.ylabel(r"$z/R_1$")
            plt.colorbar(im,label=r"radial velocity (inertial)")
            
            plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
            
            plt.savefig(output_dir+"vr_vert_"+str(i)+".png",bbox_inches='tight',dpi=100)
            plt.close()

        if 'divV' in plot_vars:

            divV = ( d['gx1v']**(-2) *  np.gradient(d['gx1v']**2 * d['vel1'],axis=1)/np.gradient(d['gx1v'],axis=1) +
                     1/(d['gx1v']*np.sin(d['gx2v'])) *  np.gradient(np.sin(d['gx2v']) *d['vel2'],axis=0)/np.gradient(d['gx2v'],axis=0) )
         
            plt.figure(figsize=(12,9) )
            im=plt.pcolormesh(x-np.linalg.norm(rcom),z,
                              divV,
                              cmap=plt.cm.PRGn,
                              rasterized=True,
                              vmin=-0.5,vmax=0.5)
        
            plt.plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
                     z2,
                     'w*',markersize=3)
            
            plt.axis('equal')
            plt.xlabel(r"$x'/R_1$")
            plt.ylabel(r"$z/R_1$")
            plt.colorbar(im,label=r"$\nabla v$")
            
            plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=1)),(-29,29),color='k',fontsize='small')
            
            plt.savefig(output_dir+"divV_vert_"+str(i)+".png",bbox_inches='tight',dpi=100)
            plt.close()
