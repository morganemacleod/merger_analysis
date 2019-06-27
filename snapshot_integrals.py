import matplotlib as mpl
mpl.use('agg')
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from glob import glob
import OrbitAnalysisUtils as ou
import argparse

### SIMULATION PARAMS ######
parser = argparse.ArgumentParser(description='Read m1,m2, input/output directories')

parser.add_argument("--m1",type=float,help="mass of particle m1",default=0.0)
parser.add_argument("--m2",type=float,help="mass of particle m2",default=0.0)

parser.add_argument("--base_dir", help="data directory (should end with / )")
#parser.add_argument("--output_dir", help="directory to save figures/output (should end with / )")
parser.add_argument("--first_file_index",help="neg number of index for first file in list eg. -100 (or 0 to start at the beginning)",default=-30,type=int)
parser.add_argument("--gamma",help="adiabatic index",default=5./3.,type=float)

args = parser.parse_args()
m1=args.m1
m2=args.m2
base_dir=args.base_dir
#output_dir=args.output_dir



filelist = sorted(glob(base_dir+"HSE.out0.00[0-9][0-9][0-9].athdf"))
#filelist = sorted(glob(base_dir+"HSE.out1.005[0-9]0.athdf"))
#filelist = filelist[args.first_file_index:]
print filelist

radii = [1,2,3,4,6,10,15,20,30]
names = ['time','sep','mass_bound','mass_unbound','mass_bound_rL1','mass_unbound_rL1','rL1','r1','r2','r3','r4','r6','r10','r15','r20','r30']
############################

orb = ou.read_trackfile(base_dir+"pm_trackfile.dat",m1=m2,m2=m2)
#print "ORB: ... ", orb.colnames

data = []

for i,myfile in enumerate(filelist):
    
    d=ou.read_data(myfile,orb,m1=m1,m2=m2,rsoft2=0.05,level=0,get_cartesian=True,get_torque=False,get_energy=True,
                   profile_file=base_dir+"hse_profile.dat",gamma=args.gamma)
    
    t = d['Time']
    sep =  np.interp(t,orb['time'],orb['sep'])
    
    data_entry = [t,sep]
    select_unbound = ((d['bern']>0) & (d['gx1v']>1.0))
    mu = np.sum(d['rho'][select_unbound]*d['dvol'][select_unbound])
    select_bound = ((d['bern']<=0) & (d['gx1v']>1.0))
    mb = np.sum(d['rho'][select_bound]*d['dvol'][select_bound])
    data_entry.append(mb)
    data_entry.append(mu)

    del select_unbound, select_bound

    ## Within L1 special case (assume m1=1)
    q = m2/1.0
    rL1 = max(1.0, sep / (1.0+np.sqrt(q)) ) # larger of between R_donor and L1 point
    select_L1 = d['gx1v']<=rL1
    mL1 = np.sum(d['rho'][select_L1]*d['dvol'][select_L1])
    select_unbound_L1 = ((d['bern']>0) & (d['gx1v']>rL1))
    muL1 = np.sum(d['rho'][select_unbound_L1]*d['dvol'][select_unbound_L1])
    select_bound_L1 = ((d['bern']<=0) & (d['gx1v']>rL1))
    mbL1 = np.sum(d['rho'][select_bound_L1]*d['dvol'][select_bound_L1])
    data_entry.append(mbL1)
    data_entry.append(muL1)
    data_entry.append(mL1)

    del select_L1, select_unbound_L1, select_bound_L1

    ## Other radii
    for my_r in radii:
        select = d['gx1v']<=my_r
        my_m  = np.sum(d['rho'][select]*d['dvol'][select])
        data_entry.append(my_m)
  
    
    data.append(data_entry)

datatable = Table(np.array(data),names=names )
ascii.write(datatable,output=base_dir+"mass_time_L1.dat")
