import numpy as np
from astropy.io import ascii
from astropy.table import Table
from glob import glob
import OrbitAnalysisUtils as ou


### SIMULATION PARAMS ######
m1 = 0.631686
m2 = 0.3

base_dir = "/Users/morganmacleod/DATA/athenaruns/pm_envelope/smr_RL_hr_lr-diode/"
output_dir = "diode_figures/"

filelist = glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf")
filelist = filelist[-41:-1]
print filelist
############################

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
#print "ORB: ... ", orb.colnames

data = []

for i,myfile in enumerate(filelist):
    d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.1,level=0,get_cartesian=True,get_torque=True, x1_max=3.5)

    t = d['Time']
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)

    theta_rot = -np.arctan2(y2-rcom[1],x2-rcom[0])

    xrotfull = (d['x']-rcom[0])*np.cos(theta_rot) - (d['y']-rcom[1])*np.sin(theta_rot)
    yrotfull = (d['x']-rcom[0])*np.sin(theta_rot) + (d['y']-rcom[1])*np.cos(theta_rot)

    xL,phiL,bin_phi = ou.get_roche_function(orb,t,M1=1,M2=m2)

    phi = bin_phi(xrotfull,yrotfull,d['z'])

    phiL1=min(phiL)

    dist2 = np.sqrt( (d['x']-x2)**2 + (d['y']-y2)**2 + (d['z']-z2)**2 )
    dist1 = d['gx1v']

    sep = np.sqrt( x2**2 + y2**2 + z2**2 )
    selectRL2 = (phi<phiL1) & (xrotfull>xL[1]) & (dist2<0.75*sep*(0.3/1.0)**(1./3.)) 
    selectRL1 = (phi<phiL1) & (xrotfull<xL[1]) & (dist1<0.75*sep*(1.0/0.3)**(1./3.)) 
    
    mRL1 = np.sum(d['rho'][selectRL1]*d['dvol'][selectRL1])
    mRL2 = np.sum(d['rho'][selectRL2]*d['dvol'][selectRL2])
    
    t2_RL1 = np.sum(d['torque_dens_2_z'][selectRL1] * d['dvol'][selectRL1])
    t2_RL2 = np.sum(d['torque_dens_2_z'][selectRL2] * d['dvol'][selectRL2])
    
    t1_RL1 = np.sum(d['torque_dens_1_z'][selectRL1] * d['dvol'][selectRL1])
    t1_RL2 = np.sum(d['torque_dens_1_z'][selectRL2] * d['dvol'][selectRL2])
    
    data.append([t,mRL1,mRL2,t1_RL1,t1_RL2,t2_RL1,t2_RL2])

datatable = Table(np.array(data),names=['time','mRL1','mRL2','t1_RL1','t1_RL2','t2_RL1','t2_RL2'] )
ascii.write(datatable,output="roche_mass_torque_time.dat")
