import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
import OrbitAnalysisUtils as ou

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
#base_dir = "/Users/morganmacleod/DATA/athenaruns/pm_envelope/smr_RL_hr_lr2/"
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/smr_RL_hr_lr2/"

output_dir = "paper_figures/"

m1 = 0.631686
m2 = 0.3
G=1

myfile = base_dir+"HSE.out1.00369.athdf"

mycm = plt.cm.RdBu_r

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





##############
# ENERGY
##############

x2slicevalue=get_midplane_theta(myfile,level=mylevel)
print "Slicing at x2=",x2slicevalue

d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.1,level=2,get_cartesian=True,get_torque=False,get_energy=True, x1_max=30,x2_min=x2slicevalue,x2_max=x2slicevalue,
               profile_file=base_dir+"hse_profile.dat")


#fig = plt.figure(figsize=(8,11))
t = d['Time']
    
rcom,vcom = ou.rcom_vcom(orb,t)
x2,y2,z2 = ou.pos_secondary(orb,t)
    
theta_rot = -np.arctan2(y2,x2)

xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)


im=plt.pcolormesh(ou.get_plot_array_midplane(xrot),
                  ou.get_plot_array_midplane(yrot),
                  ou.get_plot_array_midplane(d['etot'][:,len(d['x2v'])/2,:] /
                                          d['rho'][:,len(d['x2v'])/2,:]) ,
                  cmap=mycm,rasterized=True,
                  vmin=-1,vmax=1)
    

plt.contour(ou.get_plot_array_midplane(xrot),
               ou.get_plot_array_midplane(yrot),
               ou.get_plot_array_midplane(d['etot'][:,len(d['x2v'])/2,:] /
                                          d['rho'][:,len(d['x2v'])/2,:]) ,
               levels=[0],colors='k',linewidths=0.5)
    
plt.plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
         (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
         'w*',markersize=3)
    
plt.xlabel(r"$x'/R_1$")
plt.ylabel(r"$y'/R_1$")


plt.colorbar(im,label='specific binding energy $[G M_1 / R_1]$')

plt.xlim(-10,10)
plt.ylim(-10,10)
plt.xticks([-10,-5,0,5,10])
plt.yticks([-10,-5,0,5,10])

plt.savefig(output_dir+"binding_energy_slice_final.pdf",dpi=300,bbox_inches='tight')
plt.clf()


# VERTICAL
dblank=ar.athdf(myfile,level=mylevel,quantities=[],subsample=True)
t=dblank['Time']
rcom,vcom = ou.rcom_vcom(orb,t)
x2,y2,z2 = ou.pos_secondary(orb,t)

theta_rot = -np.arctan2(y2,x2)
x3slicevalue = dblank['x3v'][np.argmin(np.abs(dblank['x3v']+theta_rot))]
print "x3=",x3slicevalue
# read the data
d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,get_energy=True,
                 x1_max=30,x3_min=x3slicevalue,x3_max=x3slicevalue,profile_file=base_dir+"hse_profile.dat")
        
        
im=plt.pcolormesh(d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                      d['z'][0,:,:],
                      d['etot'][0,:,:]/d['rho'][0,:,:],
                      cmap=mycm,
                      vmin=-1,vmax=1,rasterized=True)

plt.contour(d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                      d['z'][0,:,:],
                      d['etot'][0,:,:]/d['rho'][0,:,:],
               levels=[0],colors='k',linewidths=0.5)


        
if(x3slicevalue<0):
    x3slicevalue += np.pi
else:
    x3slicevalue -= np.pi
print "x3=",x3slicevalue
# read the data
d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,get_energy=True,
                 x1_max=30,x3_min=x3slicevalue,x3_max=x3slicevalue,profile_file=base_dir+"hse_profile.dat")

im=plt.pcolormesh(-d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                      d['z'][0,:,:],
                      d['etot'][0,:,:]/d['rho'][0,:,:],
                      cmap=mycm,
                      vmin=-1,vmax=1,rasterized=True)

plt.contour(-d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                      d['z'][0,:,:],
                      d['etot'][0,:,:]/d['rho'][0,:,:],
               levels=[0],colors='k',linewidths=0.5)

plt.plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
             z2,
             'w*',markersize=3)


plt.colorbar(im,label='specific binding energy $[G M_1 / R_1]$')


plt.xlim(-10,10)
plt.ylim(-10,10)
plt.xticks([-10,-5,0,5,10])
plt.yticks([-10,-5,0,5,10])

plt.xlabel(r"$x'/R_1$")
plt.ylabel(r"$z/R_1$")

plt.savefig(output_dir+"binding_energy_vert_slice_final.pdf",dpi=300,bbox_inches='tight')
plt.clf()

#############################
# DENSITY SLICE
#############################
x2slicevalue=get_midplane_theta(myfile,level=mylevel)
print "Slicing at x2=",x2slicevalue

d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.1,level=2,get_cartesian=True,get_torque=False,get_energy=True, x1_max=30,x2_min=x2slicevalue,x2_max=x2slicevalue,
               profile_file=base_dir+"hse_profile.dat")

im=plt.pcolormesh(ou.get_plot_array_midplane(xrot),
                  ou.get_plot_array_midplane(yrot),
                  ou.get_plot_array_midplane( np.log10( d['rho'][:,len(d['x2v'])/2,:] ) )  ,
                  cmap=plt.cm.bone_r,rasterized=True,
                  vmin=-8,vmax=0)
    

#plt.contour(ou.get_plot_array_midplane(xrot),
#               ou.get_plot_array_midplane(yrot),
#               ou.get_plot_array_midplane(d['etot'][:,len(d['x2v'])/2,:] /
#                                          d['rho'][:,len(d['x2v'])/2,:]) ,
#               levels=[0],colors='k',linewidths=0.5)
    
plt.plot((x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot),
         (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot),
         'w*',markersize=3)
    
plt.xlabel(r"$x'/R_1$")
plt.ylabel(r"$y'/R_1$")


plt.colorbar(im,label=r'$\log \ \rho \ \ [M_1 / R_1^3]$')

plt.xlim(-25,25)
plt.ylim(-25,25)
#plt.xticks([-10,-5,0,5,10])
#plt.yticks([-10,-5,0,5,10])
plt.savefig(output_dir+"density_slice_final.pdf",dpi=300,bbox_inches='tight')
plt.clf()


# VERTICAL
dblank=ar.athdf(myfile,level=mylevel,quantities=[],subsample=True)
t=dblank['Time']
rcom,vcom = ou.rcom_vcom(orb,t)
x2,y2,z2 = ou.pos_secondary(orb,t)

theta_rot = -np.arctan2(y2,x2)
x3slicevalue = dblank['x3v'][np.argmin(np.abs(dblank['x3v']+theta_rot))]
print "x3=",x3slicevalue
# read the data
d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                 x1_max=30,x3_min=x3slicevalue,x3_max=x3slicevalue)
        
        
im=plt.pcolormesh(d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                      d['z'][0,:,:],
                      np.log10(d['rho'][0,:,:]),
                      cmap=plt.cm.bone_r,
                      vmin=-8,vmax=0,rasterized=True)
        
if(x3slicevalue<0):
    x3slicevalue += np.pi
else:
    x3slicevalue -= np.pi
print "x3=",x3slicevalue
# read the data
d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.1,level=mylevel,get_cartesian=True,get_torque=False,
                 x1_max=30,x3_min=x3slicevalue,x3_max=x3slicevalue)

im=plt.pcolormesh(-d['gx1v'][0,:,:]*np.sin(d['gx2v'][0,:,:])-np.linalg.norm(rcom),
                      d['z'][0,:,:],
                      np.log10(d['rho'][0,:,:]),
                      cmap=plt.cm.bone_r,
                      vmin=-8,vmax=0,rasterized=True)

plt.plot(np.sqrt(x2**2 + y2**2)-np.linalg.norm(rcom),
             z2,
             'w*',markersize=3)





plt.colorbar(im,label=r'$\log \ \rho \ \ [M_1 / R_1^3]$')
plt.xlabel(r"$x'/R_1$")
plt.ylabel(r"$z/R_1$")
plt.xlim(-25,25)
plt.ylim(-25,25)
plt.savefig(output_dir+"density_vert_slice_final.pdf",dpi=300,bbox_inches='tight')
plt.clf()
