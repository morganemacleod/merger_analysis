import numpy as np
from astropy.io import ascii
from astropy.table import Table
from glob import glob
import OrbitAnalysisUtils as ou


### SIMULATION PARAMS ######
m1 = 0.631686
m2 = 0.3

base_dir = "/Users/morganmacleod/DATA/athenaruns/pm_envelope/smr_RL_hr_lr/"

filelist = glob(base_dir+"HSE.out1.000[3-9][0-9].athdf")
#filelist = glob(base_dir+"HSE.out1.00039.athdf")
print filelist

radii = [1,2,3,4,6,10,15,20,30]
names = ['time','r1','r2','r3','r4','r6','r10','r15','r20','r30']
############################

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
#print "ORB: ... ", orb.colnames

data = []

for i,myfile in enumerate(filelist):
    
    d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.1,level=0,get_cartesian=True,get_torque=False)
    
    t = d['Time']
    
    data_entry = [t]
    for my_r in radii:
        rindex = np.argmin( np.abs(d['x1v']-my_r))
        my_m  = np.sum(d['rho'][rindex,:,:]*d['dA'][rindex,:,:]*d['vel1'][rindex,:,:])
        data_entry.append(my_m)
  
    
    data.append(data_entry)

datatable = Table(np.array(data),names=names )
ascii.write(datatable,output="mdot_radii_time.dat")

