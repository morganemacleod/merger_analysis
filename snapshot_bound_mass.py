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

#filelist = glob(base_dir+"HSE.out1.002[0-9][0-9].athdf")
filelist = glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf")
filelist = filelist[-41:-1]
print filelist

radii = [1,2,3,4,6,10,15,20,30]
names = ['time','mass_bound','mass_unbound']
############################

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
#print "ORB: ... ", orb.colnames

data = []

for i,myfile in enumerate(filelist):
    
    d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.1,level=0,get_cartesian=True,get_torque=False,get_energy=True,
                  profile_file=base_dir+"hse_profile.dat")
    
    t = d['Time']
    
    data_entry = [t]
    select_unbound = d['etot']>0
    mu = np.sum(d['rho'][select_unbound]*d['dvol'][select_unbound])
    select_bound = d['etot']<=0
    mb = np.sum(d['rho'][select_bound]*d['dvol'][select_bound])
    data_entry.append(mb)
    data_entry.append(mu)
  
    
    data.append(data_entry)

datatable = Table(np.array(data),names=names )
ascii.write(datatable,output=output_dir+"mass_bound_time.dat")
