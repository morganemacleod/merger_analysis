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

parser.add_argument("m1",type=float,help="mass of particle m1")
parser.add_argument("m2",type=float,help="mass of particle m2")

parser.add_argument("--base_dir", help="data directory (should end with / )")
parser.add_argument("--output_dir", help="directory to save figures/output (should end with / )")
parser.add_argument("--first_file_index",help="neg number of index for first file in list eg. -100",default=-30,type=int)
parser.add_argument("--last_file_index",help="neg number of index for first file in list eg. -1",default=-1,type=int)


args = parser.parse_args()
m1=args.m1
m2=args.m2
base_dir=args.base_dir
output_dir=args.output_dir


filelist = sorted(glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf"))
filelist = filelist[args.first_file_index:args.last_file_index]
print filelist

names = ['time','mdot']
############################

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
#print "ORB: ... ", orb.colnames


# GET THE SLICE VALUE
slice_radius = 1.0
dblank = ou.ar.athdf(filelist[0],level=2,quantities=[],subsample=True)
x1slice = dblank['x1v'][np.argmin( np.abs(dblank['x1v']-slice_radius) )]

# LOOP THROUGH FILES
data = []

for i,myfile in enumerate(filelist):

    d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.05,level=0,get_cartesian=True,get_torque=False,get_energy=False,
                   profile_file=base_dir+"hse_profile.dat",x1_min=x1slice,x1_max=x1slice)
    
    t = d['Time']

    md = np.sum(d['dA']*d['rho']*d['vel1'])
    data_entry = [t,md]
    
    data.append(data_entry)

datatable = Table(np.array(data),names=names )
ascii.write(datatable,output=output_dir+"mdot_time.dat")
