import numpy as np
from astropy.io import ascii
from astropy.table import Table
from glob import glob
import OrbitAnalysisUtils as ou


### SIMULATION PARAMS ######
m1 = 0.410103
m2 = 0.3

base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/syncRL/fcorot_series/a206_res24_fc10/"

filelist = glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf")
filelist = sorted(filelist)[-100:]
#filelist = glob(base_dir+"HSE.out1.00539.athdf")
print filelist

radii = [ 0.1,  0.2,  0.3,  0.4,  0.5,  0.6 ]
names = ['time','r01','r02','r03','r04','r05','r06']
############################

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
#print "ORB: ... ", orb.colnames

data = []

for i,myfile in enumerate(filelist):
    
    d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.1,level=0,get_cartesian=True,get_torque=True,x1_max=3.25)
    
    t = d['Time']
    x2,y2,z2 = ou.pos_secondary(orb,t)
    
    data_entry = [t]
    for my_r in radii:
        d2 = np.sqrt((d['x']-x2)**2 + (d['y']-y2)**2 + (d['z']-z2)**2 )
        select = d2<=my_r
        my_t  = np.sum(d['torque_dens_2_z'][select]*d['dvol'][select])
        data_entry.append(my_t)
  
    
    data.append(data_entry)

datatable = Table(np.array(data),names=names )
ascii.write(datatable,output="torque2_radii_time.dat")

