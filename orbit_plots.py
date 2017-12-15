# normal stuff
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
import OrbitAnalysisUtils as ou
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


def a_RL(q):
    """Eggelton formula, q=M2/M1 (opposite of eggelton deff)"""
    return (0.6*q**(-2./3.) + np.log(1+q**(-1./3.)))/(0.49*q**(-2./3.))

def p_orb(q,a):
    return 2*np.pi*np.sqrt(a**3/(1+q))

####### SIMULATION PARAMS  ############
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

#######################################

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
print "ORB: ... ", orb.colnames

hst = ascii.read(base_dir+"HSE.hst",
                 names=['time','dt','mass','1-mom','2-mom','3-mom','1-KE','2-KE','3-KE','tot-E','mxOmegaEnv','mEnv'])
print "\nHSE: ...", hst.colnames

mg = hst['mass'][0]
print mg+m1

t1=ou.get_t1(orb)
print "t1=",t1


####### 
# ORBITAL PATH FIGS
#######
plt.figure(figsize=(5,5))
phi = np.linspace(0,2*np.pi,100)
plt.plot(orb['x'],orb['y'],'k-',label="$m_2$")
plt.plot(np.cos(phi),np.sin(phi),'--',color='grey')
plt.plot(0.3*np.cos(phi),0.3*np.sin(phi),'-',color='grey')
plt.plot(0,0,'bx',label="$m_1$")
plt.axis('equal')
plt.title('Simulation frame')
plt.xlabel(r'$x/R_1$')
plt.ylabel(r'$y/R_1$')
plt.xlim(-2.5,2.5)
plt.ylim(-2.5,2.5)
plt.legend(loc=0)
plt.savefig(output_dir+"orbital_path_sim.pdf",bbox_inches='tight')
plt.clf()

plt.figure(figsize=(5,5))
phi = np.linspace(0,2*np.pi,100)
plt.plot(orb['x']-orb['xcom'],orb['y']-orb['ycom'],'k-',label="$m_2$")
plt.plot(-orb['xcom'],-orb['ycom'],'b-',label="$m_1$")
plt.title('Inertial frame')
plt.xlabel(r'$x/R_1$')
plt.ylabel(r'$y/R_1$')
plt.axis('equal')
plt.xlim(-2.5,2.5)
plt.ylim(-2.5,2.5)
plt.legend(loc=0)
plt.savefig(output_dir+"orbital_path_in.pdf",bbox_inches='tight')
plt.clf()

plt.figure(figsize=(9,4))
plt.plot(orb['time']-t1,orb['sep'],'k-')
plt.axhline(1,ls='--',color='grey')
plt.axhline(0.3,ls='-',color='grey')
plt.axhline(a_RL(0.3),ls='-',color='SkyBlue')
plt.ylabel(r'separation $[R_1]$')
plt.xlabel(r'$t-t_1 \ [(R_1^3/GM_1)^{1/2}]$')
plt.xlim(orb['time'][0]-t1,10)
plt.savefig(output_dir+"separation_time.pdf",bbox_inches='tight')
plt.clf()



###### ANGULAR MOMENTUM 
plt.figure()
plt.plot(orb['time']-t1, orb['lpz'], label=r'particle $L_z$',lw=2,color='CadetBlue'  )
plt.plot(orb['time']-t1, orb['lgz'], label=r'gas $L_z$',ls='--',lw=2,color='RosyBrown'  )
plt.plot(orb['time']-t1, orb['lgoz'],label=r'gas off grid $L_z$',ls=':',lw=2,color='RosyBrown' )

plt.plot(orb['time']-t1, orb['ltz'],label=r'total $L_z$',lw=2,color='k' )
plt.legend(loc=0)
plt.xlabel(r'$t-t_1 \ [(R_1^3/GM_1)^{1/2}]$')
plt.ylabel(r'$\hat z$ angular momenta $\left[(GM_1 R_1)^{1/2} \right]$')
plt.xlim(orb['time'][0]-t1,10)
plt.savefig("paper_figures/ang_mom_time.pdf",bbox_inches='tight')
plt.clf()

plt.plot(orb['sep'], orb['ltz']/orb['ltz'][0],label=r'total $L_z$',lw=2,color='k' )
plt.legend(loc=0)
plt.xlabel(r'separation $[R_1]$')
#plt.ylabel(r'$\hat z$ angular momenta $\left[(GM_1 R_1)^{1/2} \right]$')
plt.savefig(output_dir+"norm_tot_ang_mom_sep.pdf",bbox_inches='tight')
plt.clf()


##### OMEGA 
plt.figure()
# an approximate version of omega -- this is strictly valid for a circular orbit but should be pretty close
Omega = orb['vmag']/orb['sep']
plt.plot(orb['sep'],Omega/Omega[0],label=r"$\Omega_{\rm orb}$",color='Black',lw=2)

Omega_env = hst['mxOmegaEnv']/hst['mEnv']
plt.plot(np.interp(hst['time'],orb['time'],orb['sep']),Omega_env/Omega[0],'--',label=r"$\langle \Omega_{\rm env} \rangle$",color='RosyBrown',lw=2)

plt.plot(orb['sep'],(orb['sep']/orb['sep'][0])**-1.5,':',lw=2,color='grey',label=r'$a^{-3/2}$'  )

#plt.xlabel(r'time $[R_1^3/GM_1]$')
plt.ylabel(r"$\Omega/\Omega_{\rm orb,0}$")
plt.xlabel(r'separation $[R_1]$')
plt.legend(loc=0)
plt.savefig(output_dir+"omega_separation.pdf",bbox_inches='tight')
plt.clf()


### Decay timescale
plt.figure()
# number of orbits for decay of angular momentum
Omega = orb['vmag']/orb['sep']
dadt = np.gradient(orb['sep'])/np.gradient(orb['time'])

#plot_arr = convolve( orb['sep']/np.abs(dadt) * Omega/(2.*np.pi),  Box1DKernel(1) )
plot_arr = orb['sep']/np.abs(dadt) * Omega/(2.*np.pi)

plt.plot(orb['sep'],plot_arr,
         'k-',label=r'$a\Omega/(2\pi |\dot a|) $' )

plt.legend(loc=0)

plt.yscale('log')
plt.xlabel(r'separation $[R_1]$')
plt.ylabel(r'$N_{\rm decay}$')
plt.ylim(0.3,1.e2)
plt.xlim(0.6,2.06)
plt.savefig(output_dir+"Ndecay_sep.pdf",bbox_inches='tight')
plt.clf()




########
# Torque plot
#######

# TOTAL GRAV TORQUE ON THE TWO COMPONENTS
orb['t2'] = np.cross(orb['r']-orb['rcom'],m2*orb['agas2'])[:,2]
orb['t1'] = np.cross(-orb['rcom'],m1*orb['agas1'])[:,2]

tmin=-30
select = orb['time']-t1>tmin

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.plot(orb[select]['time']-t1,orb[select]['t1'],'--',label=r'$1$',lw=2)
plt.plot(orb[select]['time']-t1,orb[select]['t2'],':',label=r'$2$',lw=2)
plt.plot(orb[select]['time']-t1,orb[select]['t1']+orb[select]['t2'],'k-',label='total')
plt.legend(loc=0,handletextpad=0.5,frameon=True)
plt.xlabel("$t-t_1 \ \ [( R_1^3 / GM_1 )^{1/2}]$")
plt.ylabel(r'$\tau_{\rm grav}$')
plt.grid()

plt.subplot(122)
plt.plot(orb[select]['sep'],orb[select]['t1'],'--',label=r'$1$',lw=2)
plt.plot(orb[select]['sep'],orb[select]['t2'],':',label=r'$2$',lw=2)
plt.plot(orb[select]['sep'],orb[select]['t1']+orb[select]['t2'],'k-',label='total')
#plt.legend(loc=0,handletextpad=0.5)
plt.xlabel(r'Separation $[R_1]$')
#plt.ylabel(r'$\tau_{\rm grav}$')
plt.yticks(visible=False)
plt.grid()

plt.subplots_adjust(wspace=0.0)
plt.savefig(output_dir+"torque_time_sep.pdf",bbox_inches='tight')
