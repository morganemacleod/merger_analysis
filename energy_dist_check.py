import numpy as np
import matplotlib.pyplot as plt
import athena_read as ar
import OrbitAnalysisUtils as ou
from Constants import Constants
import argparse
#import seaborn as sns
#import deepdish as dd
#from astropy.table import Table
#from glob import glob
#from mpl_toolkits.axes_grid1 import ImageGrid
#from tqdm.auto import tqdm

c=Constants()

#%matplotlib inline

plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.2
plt.rcParams['legend.labelspacing'] = 0.2
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 16


parser = argparse.ArgumentParser(description='Read filename to process')

parser.add_argument("--filename", help="full 3d file (athdf)")

args = parser.parse_args()



base_dir = "./"

orb = ou.read_trackfile(base_dir+"pm_trackfile.dat")

rin = 2.45e11
plt.figure(figsize=(8,8))
plt.subplot(221)
plt.plot(orb['time'],orb['sep'])
plt.xlabel('t')
plt.ylabel('sep')
plt.axhline(rin,color='grey')
plt.semilogy()

plt.subplot(222)
plt.plot(orb['x'],orb['y'])
pp = np.linspace(0,2*np.pi,1000)
plt.plot( rin*np.cos(pp), rin*np.sin(pp),'r-' )
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')

plt.subplot(223)
plt.plot(orb['time'],orb['m1'])
plt.plot(orb['time'],orb['m2'])
plt.plot(orb['time'],orb['m1']+orb['m2'])
plt.xlabel('t')
plt.ylabel('mass')
plt.axhline(rin,color='grey')

plt.subplot(224)
plt.plot(orb['time'],orb['dt'])
plt.xlabel('t')
plt.ylabel('dt')
plt.semilogy()

plt.subplots_adjust(hspace=0.25,wspace=0.25)
plt.savefig("orb_checks.png",bbox_inches='tight')





### 3D READ
base_dir = "./"

orb = ou.read_trackfile(base_dir+"pm_trackfile.dat")
myfile = base_dir + args.filename 
lim= 40
thind = 0
mylevel = 0
r0thresh = 0.
velvecs = True
skip = 4
x2slicevalue=ou.get_midplane_theta(myfile,level=mylevel)
d = ou.read_data(myfile,orb,gamma=1.6,level=mylevel,rsoft2=3.25e10,
                   get_energy=True,profile_file=base_dir+'hse_profile.dat')
x2,y2,z2 = ou.pos_secondary(orb,d['Time'])
rcom,vcom = ou.rcom_vcom(orb,d['Time'])
print("t=",d['Time'])
print('sep/rstar = ',np.sqrt(x2**2 + y2**2 + z2**2)/2.45e12)
print('total mass = ',np.sum(d['rho']*d['dvol']*d['r0'])/c.msun)
print('unbound mass =',np.sum( (d['rho']*d['dvol']*d['r0'])[d['etot']>0] ) / c.msun)
print('total mass + m1 ',(np.sum(d['rho']*d['dvol']*d['r0'])+ orb['m1'][0])/c.msun )

plt.figure()
plt.hist( d['vel1'].flatten(),weights=(d['rho']*d['dvol']).flatten(), histtype='step', lw=2,bins=30 )
plt.hist( d['vel1'].flatten(),weights=(d['rho']*d['dvol']*d['r0']).flatten(), histtype='step', lw=2,bins=30 )
plt.semilogy()
vesc = np.sqrt(2*c.G*8e33/2.45e12)
plt.axvline(0)
plt.axvline(vesc)
xp = 1e7*np.linspace(0,9)
plt.plot(xp,0.01*c.msun*np.exp(-(xp/vesc)**3 ) )
#plt.plot(xp,0.1*c.msun*np.exp(-(6*xp/vesc) ) )
plt.ylim(1e25,)
plt.xlabel('radial velocity') 
plt.ylabel('mass')
plt.savefig("vel1_dist.png",bbox_inches='tight')

plt.figure()
plt.hist((d['etot']/d['rho']).flatten(),weights=(d['rho']*d['dvol']).flatten(), histtype='step', lw=2,bins=30 )
plt.hist((d['ek']/d['rho']).flatten(),weights=(d['rho']*d['dvol']).flatten(), histtype='step', lw=2,bins=30 )
xp = 3e14*np.linspace(0,1.5)
plt.plot(xp,0.01*c.msun*np.exp(-(xp/4e13) ) )
plt.semilogy()
plt.xlabel('energy') 
plt.ylabel('mass')
plt.savefig("energy_dist.png",bbox_inches='tight')
