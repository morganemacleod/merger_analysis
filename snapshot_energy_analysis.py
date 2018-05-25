import numpy as np
from astropy.io import ascii
from astropy.table import Table
from glob import glob
import OrbitAnalysisUtils as ou

### SIMULATION PARAMS ######
m1 = 0.410103
m2 = 0.3

base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/syncRL/ejecta_v1/a206_res24_fc10/"
output_dir = ""

#filelist = glob(base_dir+"HSE.out1.002[0-9][0-9].athdf")
#filelist = glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf")
filelist = glob(base_dir+"HSE.out1.005[0-9][0,5].athdf")
#filelist = filelist[-41:-1]
print filelist

names = ['time','EK1','EK2','EPp','EPg','EPpg','EKg','EIg','Egas','Etot','EPg_sg','Etot_sg']
############################

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
#print "ORB: ... ", orb.colnames




def epotSG(d):
   """ return the gas self-gravitational potential """
   dmr = np.sum(d['rho']*d['dvol'],axis=(0,1))
   mr = np.cumsum(dmr)
   return np.sum(-(mr)*dmr/d['x1v'])




data = []

for i,myfile in enumerate(filelist):
    
    d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.05,level=0,get_cartesian=True,get_torque=False,get_energy=True,
                   profile_file=base_dir+"hse_profile.dat",x1_max=30)
    
    t = d['Time']
    
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)

    sep = np.sqrt(x2**2 + y2**2 + z2**2)
    vx2 = np.interp(t,orb['time'],orb['vx']) -vcom[0]
    vy2 = np.interp(t,orb['time'],orb['vy']) -vcom[1]
    vz2 = np.interp(t,orb['time'],orb['vz']) -vcom[2]

    EK1 =  0.5*m1*(vcom[0]**2 + vcom[1]**2 + vcom[2]**2)
    EK2 = 0.5*m2*(vx2**2 + vy2**2 + vz2**2)
    EPp = -m1*m2/sep
    

    EPg = np.sum(d['epotg']*d['dvol'])
    EPpg = np.sum(d['epotp']*d['dvol'])
    EKg = np.sum(d['ek']*d['dvol'])
    EIg = np.sum(d['ei']*d['dvol'])
    Egas = np.sum(d['etot']*d['dvol'])

    EPg_sg = epotSG(d)
    
    

    data_entry = [t]
    data_entry.append(EK1)
    data_entry.append(EK2)
    data_entry.append(EPp)
    data_entry.append(EPg)
    data_entry.append(EPpg)
    data_entry.append(EKg)
    data_entry.append(EIg)
    data_entry.append(Egas)
    data_entry.append(Egas+EK1+EK2+EPp)
    data_entry.append(EPg_sg)
    data_entry.append(EPg_sg + EPpg + EKg + EIg + EK1 + EK2 + EPp)
    
    data.append(data_entry)

datatable = Table(np.array(data),names=names )
ascii.write(datatable,output=output_dir+"energies_time.dat")
