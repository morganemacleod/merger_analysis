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



##########################
# SIMULATION PARAMS
##########################
m1 = 0.631686
m2 = 0.3

base_dir = "/Users/morganmacleod/DATA/athenaruns/pm_envelope/smr_RL_hr_lr2/"

output_dir = "paper_figures/"

orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")

t1=ou.get_t1(orb)





###
# MASS IN ROCHE LOBES
###
mt = ascii.read("roche_mass_torque_time.dat")
mt['sep'] = np.interp(mt['time'],orb['time'],orb['sep'])
mt.colnames

plt.figure(figsize=(8,5))

tmin=-30
select = mt['time']-t1>tmin

plt.subplot(121)
plt.plot(mt[select]['time']-t1,mt[select]['mRL1'],label=r'$m({\rm RL_1})$')
plt.plot(mt[select]['time']-t1,mt[select]['mRL2'],'--',label=r'$m({\rm RL_2})$')
plt.yscale('log')
plt.grid()
plt.ylim(1.e-3,1)
plt.xlim(tmin,)
plt.legend(loc=0,frameon=True)
plt.ylabel("Enclosed Gas Mass $[M_1]$")
plt.xlabel("$t-t_1 \ \ [( R_1^3 / GM_1 )^{1/2}]$")


plt.subplot(122)
plt.plot(mt[select]['sep'],mt[select]['mRL1'],label=r'$m({\rm RL_1})$')
plt.plot(mt[select]['sep'],mt[select]['mRL2'],'--',label=r'$m({\rm RL_2})$')
plt.yscale('log')
plt.yticks(visible=False)
plt.grid()
plt.ylim(1.e-3,1)
plt.xlabel("Separation $[R_1]$")

plt.subplots_adjust(wspace=0.0)
plt.savefig(output_dir+"mass_roche_time_sep.pdf",bbox_inches='tight')





###
# MASS IN RADII
###
mr = ascii.read("mass_radii_time.dat")
mr['sep'] = np.interp(mr['time'],orb['time'],orb['sep'])


plt.clf()
plt.figure(figsize=(8,5))
plt.subplot(121)
mycm = plt.cm.viridis
labels=[r"$r<1R_1$",r"$r<2R_1$",r"$r<3R_1$",r"$r<4R_1$",r"$r<6R_1$",r"$r<10R_1$",r"$r<15R_1$",r"$r<20R_1$",r"$r<30R_1$",]

tmin=-30
select = mr['time']-t1>tmin

cols = mr.colnames[1:-1]
for i,colname in enumerate(cols):
    plt.plot(mr[select]['time']-t1,mr[select][colname],
            color=mycm(float(i)/(len(cols)-1) ),
            label=labels[i],lw=2)


plt.xlabel("$t-t_1 \ \ [( R_1^3 / GM_1 )^{1/2}]$")
plt.ylabel("Enclosed Gas Mass $[M_1]$")
plt.legend(loc=0,frameon=True)
plt.grid()

plt.subplot(122)
mycm = plt.cm.viridis
labels=[r"$r<1R_1$",r"$r<2R_1$",r"$r<3R_1$",r"$r<4R_1$",r"$r<6R_1$",r"$r<10R_1$",r"$r<15R_1$",r"$r<20R_1$",r"$r<30R_1$",]

cols = mr.colnames[1:-1]
for i,colname in enumerate(cols):
    plt.plot(mr[select]['sep'],mr[select][colname],
            color=mycm(float(i)/(len(cols)-1) ),
            label=labels[i],lw=2)

plt.xlabel("Separation $[R_1]$")
plt.yticks(visible=False)
#plt.legend(loc=0)
plt.grid()


plt.subplots_adjust(wspace=0.0)
plt.savefig(output_dir+"mass_radii_time_sep.pdf",bbox_inches='tight')




###
# Unbound mass
###

mb = ascii.read("mass_bound_time.dat")
mb['sep'] = np.interp(mb['time'],orb['time'],orb['sep'])


tmin=-30
select = mb['time']-t1>tmin

plt.figure(figsize=(8,5))
plt.subplot(121)
plt.plot(mb[select]['time']-t1,mb[select]['mass_unbound'],'kx-')
plt.xlabel("$t-t_1 \ \ [( R_1^3 / GM_1 )^{1/2}]$")
plt.ylabel("Unbound Mass $[M_1]$")

plt.subplot(122)
plt.plot(mb[select]['sep'],mb[select]['mass_unbound'],'kx-')
plt.xlabel("Separation $[R_1]$")
plt.yticks(visible=False)

plt.subplots_adjust(wspace=0.0)
plt.savefig("paper_figures/mass_bound_time_sep.pdf",bbox_inches='tight')



###
# Torque spatial decomposition
###

mt = ascii.read("roche_mass_torque_time.dat")
mt['sep'] = np.interp(mt['time'],orb['time'],orb['sep'])

# TOTAL GRAV TORQUE ON THE TWO COMPONENTS
orb['t2'] = np.cross(orb['r']-orb['rcom'],m2*orb['agas2'])[:,2]
orb['t1'] = np.cross(-orb['rcom'],m1*orb['agas1'])[:,2]

tmin=-1e10
selectorb = orb['time']-t1>tmin
select = mt['time']-t1>tmin

ymin=-0.04
ymax = 0.008


plt.figure(figsize=(12,5))

plt.subplot(131)
plt.plot(orb[selectorb]['sep'],orb[selectorb]['t1'],'C0-',lw=3,
         label=r'$\tau_{\rm grav,1}$')
plt.plot(mt[select]['sep'],mt[select]['t1_RL1'],'C0:',lw=2,
        label=r'$\tau_{\rm grav,1}({\rm RL}_1)$')
plt.plot(mt[select]['sep'],mt[select]['t1_RL2'],'C0--',lw=2,
        label=r'$\tau_{\rm grav,1}({\rm RL}_2)$')
plt.legend(loc=0,fontsize=18,frameon=True)
plt.grid()
plt.ylim(ymin,ymax)
plt.xlabel("Separation $[R_1]$")
plt.ylabel(r'$\tau_{\rm grav}$',fontsize=24)
plt.title("$m_1$")


plt.subplot(132)
plt.plot(orb[selectorb]['sep'],orb[selectorb]['t2'],'C1-',lw=3,
        label=r'$\tau_{\rm grav,2}$')
plt.plot(mt[select]['sep'],mt[select]['t2_RL1'],'C1:',lw=2,
        label=r'$\tau_{\rm grav,2}({\rm RL}_1)$')
plt.plot(mt[select]['sep'],mt[select]['t2_RL2'],'C1--',lw=2,
        label=r'$\tau_{\rm grav,2}({\rm RL}_2)$')
plt.legend(loc=0,fontsize=18,frameon=True)
plt.grid()
plt.xlabel("Separation $[R_1]$")
plt.yticks(visible=False)
plt.ylim(ymin,ymax)
plt.title("$m_2$")



plt.subplot(133)
plt.plot(orb[selectorb]['sep'],orb[selectorb]['t1']+orb[selectorb]['t2'],'k-',lw=3,
         label=r'$\tau_{\rm grav}$')
plt.plot(mt[select]['sep'],mt[select]['t1_RL1']+mt[select]['t2_RL1'],'k:',
         label=r'$\tau_{\rm grav}({\rm RL}_1)$',lw=2)
plt.plot(mt[select]['sep'],mt[select]['t1_RL2']+mt[select]['t2_RL2'],'k--',
         label=r'$\tau_{\rm grav}({\rm RL}_2)$',lw=2)
plt.legend(loc=0,fontsize=18,frameon=True)
plt.ylim(ymin,ymax)
plt.grid()
plt.yticks(visible=False)
plt.xlabel("Separation $[R_1]$")
plt.title("total")

plt.subplots_adjust(hspace=0,wspace=0)
plt.savefig(output_dir+"torque_roche_sep.pdf",bbox_inches='tight')


