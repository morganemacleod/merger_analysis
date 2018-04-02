import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
import OrbitAnalysisUtils as ou
from scipy.integrate import odeint

#%matplotlib inline

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

from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel

def a_RL(q):
    """Eggelton formula, q=M2/M1 (opposite of eggelton deff) 
    this is: a/rL where a is separation and rL is the Roche Lobe radius"""
    return (0.6*q**(-2./3.) + np.log(1+q**(-1./3.)))/(0.49*q**(-2./3.))

def p_orb(Mtot,a):
    """Orbital period"""
    return 2*np.pi*np.sqrt(a**3/Mtot)

def mdot_donor(a,m1,m2,r1):
    """Pols eqn 7.5"""
    rL = a_RL(m2/m1)**-1 * a
    A = 1.0
    mdot = -A*m1/p_orb(m1+m2,a)*((r1-rL)/r1)**3
    return np.select([mdot<0,mdot>=0],[mdot,0])

def derivs(vec,t,mode):
    """Pols Ch 7, Lajoie&Sills 2011"""
    # unpack    
    m1 = vec[0]
    a = vec[1]
    
    # constants
    r1 = 1.04
    m2 = 0.3
    if mode=='donor':
        gamma = m2/m1  # pols
    if mode=='accretor':
        gamma = 2*m1/m2  # pols
    if mode=='l2':
        gamma = (m1+m2)**2/(m1*m2) * 1.2**2  # pribula 
    
    # derivatives
    dm1dt = mdot_donor(a,m1,m2,r1)
    dadt = -2.0*a*dm1dt/m1*(1.0-(gamma+0.5)*m1/(m1+m2))
    
    return np.array([dm1dt,dadt])


m1_initial = 1.0
a_initial = 2.05
ic = np.array([m1_initial,a_initial])
t = np.linspace(0,1000,100000)

sol_acc=odeint(derivs,ic,t,args=('l2',))
sol_acc=Table(data=sol_acc,names=['m1','a'])
sol_acc['t']=t

select = sol_acc['m1']>0.5
#plt.plot(sol_acc[select]['t'],sol_acc[select]['a'])
#plt.plot(sol_acc[select]['t'],sol_acc[select]['m1'])


orb = ou.read_trackfile(0.410103,0.3,"/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/pm_trackfile_a206_rin03_r005_n24_fc10.dat")
t1=ou.get_t1(orb)

md = ascii.read("mdot_time.dat")
md['sep'] = np.interp(md['time'],orb['time'],orb['sep'])

#plt.plot(md['sep'],md['mdot'])

mdsmooth = convolve(md['mdot'],Box1DKernel(3),boundary='extend')
#plt.plot(md['sep'],mdsmooth)

select = md['sep']>1.6
mdsmooth[select] = convolve(mdsmooth[select],Box1DKernel(10),boundary='extend')
#plt.plot(md['sep'],mdsmooth)

select = md['sep']>1.8
mdsmooth[select] = convolve(mdsmooth[select],Box1DKernel(20),boundary='extend')
#plt.plot(md['sep'],mdsmooth)

#plt.yscale('log')
#plt.show()

dlgdt = np.gradient(orb['lgz'])/np.gradient(orb['time'])
#plt.plot(orb['sep'],dlgdt)

ldsmooth = convolve(dlgdt,Box1DKernel(100),boundary='extend')
#plt.plot(orb['sep'],ldsmooth)

select = orb['sep']>1.5
ldsmooth[select] = convolve(ldsmooth[select],Box1DKernel(600),boundary='extend')
#plt.plot(orb['sep'],ldsmooth)

select = orb['sep']>1.8
ldsmooth[select] = convolve(ldsmooth[select],Box1DKernel(2000),boundary='extend')
#plt.plot(orb['sep'],ldsmooth)
#plt.yscale('log')
#plt.show()


dmdtsmooth = np.interp(orb['time'],md['time'],mdsmooth)

select = ((orb['sep']<2.01) & (orb['sep']>0.6))
#plt.plot(orb['sep'][select],ldsmooth[select]/dmdtsmooth[select])

#plt.ylim(0,3)

xp = np.linspace(1,2.02,1001)

yp = np.interp(xp,np.flipud(orb['sep'][select]),np.flipud(ldsmooth[select]/dmdtsmooth[select]) )

#plt.plot(xp,yp,'-')

yps = convolve(yp,Box1DKernel(100),boundary='extend')
#plt.plot(xp,yps,'--' )


# MASS TRANSFER
plt.figure(figsize=(7,3.5))

select = (sol_acc['m1']>0.5) & (sol_acc['a']>0.6)

plt.subplot(121)
plt.plot(md['time']-t1,md['mdot'],'-',color='grey',lw=0.5 ,label='')
plt.plot(sol_acc[select]['t']-sol_acc[select]['t'][-1],
         -np.gradient(sol_acc[select]['m1'])/np.gradient(sol_acc[select]['t']),'C1--',
        label='analytic')
plt.plot(md['time']-t1,mdsmooth,'C0-',lw=2,label='simulation' )
plt.yscale('log')
plt.xlim(-250,)
plt.ylim(1.e-6,1.e-1)
plt.grid()
plt.ylabel(r"$|\dot M_{\rm donor}|$")
plt.xlabel("time, $t-t_1$")
plt.legend(loc=0,frameon=True,fontsize=16)


plt.subplot(122)
plt.plot(md['sep'],md['mdot'],'-',color='grey',lw=0.5 )
plt.plot(sol_acc[select]['a'],
         -np.gradient(sol_acc[select]['m1'])/np.gradient(sol_acc[select]['t']),'C1--',
         label='analytic' )
plt.plot(md['sep'],mdsmooth,'C0-',lw=2,label='simulation')
plt.xlabel(r'separation')
plt.yticks(visible=False)
plt.ylim(1.e-6,1.e-1)
plt.grid()
plt.yscale('log')


plt.subplots_adjust(wspace=0.0)
plt.savefig("paper_figures/mdot_time_sep.pdf",bbox_inches='tight')
plt.close()


# SPECIFIC MOMENTUM, GAMMA
plt.figure(figsize=(6,6))

m1 = 1.0 - np.cumsum(md['mdot']*np.gradient(md['time']))
print m1
m2 = 0.3
Mtot = m1+m2
mu = m1*m2/Mtot
Lbin = mu*np.sqrt(Mtot * md['sep'])

MtotoLbin = np.interp(xp,np.flipud(md['sep']),np.flipud(Mtot/Lbin) )

gamma = yp*MtotoLbin
gammas = yps*MtotoLbin

plt.subplot(211)
plt.plot(xp,yp,lw=0.5,color='grey')
plt.plot(xp,yps,lw=1.5)
plt.plot(md['sep'],m1/m2*(Lbin/Mtot),linestyle='--',color='grey')
plt.plot(md['sep'],(m1+m2)**2/(m1*m2) * 1.2**2 *(Lbin/Mtot),linestyle='-.',color='grey')
plt.xticks(visible=False)
plt.xlim(1,)
plt.ylim(0.25,2.75)
plt.ylabel(r"$h_{\rm loss}$")
plt.annotate(r"$h_{\rm accretor}$",(1.4,0.75))
plt.annotate(r"$h_{\rm L_2}$",(1.6,2.2))


plt.subplot(212)
plt.plot(xp,gamma,lw=0.5,color='grey')
plt.plot(xp,gammas,lw=1.5)
plt.plot(md['sep'],m1/m2,linestyle='--',color='grey')
plt.plot(md['sep'],(m1+m2)**2/(m1*m2) * 1.2**2,linestyle='-.',color='grey')
plt.xlabel(r'separation')
plt.xlim(1,)
plt.ylim(2.5,9.5)
plt.ylabel(r"$\gamma_{\rm loss}$")
plt.annotate(r"$\gamma_{\rm accretor}$",(1.2,3.4))
plt.annotate(r"$\gamma_{\rm L_2}$",(1.6,8.25))

plt.subplots_adjust(hspace=0.0)


plt.savefig("paper_figures/specific_ang_momentum_sep.pdf",bbox_inches='tight')
plt.close()


plt.figure(figsize=(6,3.))
t1 = ou.get_t1(orb)
t1_analytic = np.interp(1.0,np.flipud(sol_acc[select]['a']),np.flipud(sol_acc[select]['t']) )

plt.axhline(1,ls='--',color='grey')
plt.axhline(0.3,ls='-',color='grey')
plt.axhline(a_RL(0.3),ls='-',color='SkyBlue')

select = (sol_acc['m1']>0.5) & (sol_acc['a']>0.6)
plt.plot(sol_acc[select]['t']-t1_analytic,
         sol_acc[select]['a'],'C1--',label='analytic',lw=2)

plt.plot(orb['time']-t1,orb['sep'],'C0-',lw=2,label='simulation')
plt.xlabel("time, $t-t_1$")
plt.ylabel(r'separation')
plt.xlim(-250,5)
#plt.grid()

plt.legend(loc=0,frameon=True)

plt.savefig("paper_figures/separation_time_analytic_comparison.pdf",bbox_inches='tight')
plt.close()



### MDOT TIME

t1 = ou.get_t1(orb)
tmin = -40
tmax = 2

plt.figure(figsize=(6,9.))

plt.subplot(311)
plt.plot(orb['time']-t1,orb['sep'],'C0-',lw=2,label='simulation')
plt.ylabel(r'separation')
plt.xlim(tmin,tmax)
plt.ylim(0.5,2.06)
plt.grid()
plt.xticks(visible=False)

plt.subplot(312)
plt.plot(md['time']-t1,md['mdot'],'C0-',lw=2 ,label='')
#plt.yscale('log')
plt.xlim(tmin,tmax)
#plt.ylim(1.e-3,1.e-1)
plt.grid()
plt.ylabel(r"$|\dot M_{\rm donor}|$")
#plt.xlabel("$t-t_1 \ \ [( R_1^3 / GM_1 )^{1/2}]$")
plt.xticks(visible=False)

plt.subplot(313)
plt.plot(md['time']-t1,np.log10(md['mdot']),'C0-',lw=2 ,label='')
#plt.plot(md['time']-t1,mdsmooth,'C1-',lw=2 ,label='')
#plt.yscale('log')
plt.xlim(tmin,tmax)
plt.ylim(-4,-1)
plt.grid()
plt.ylabel(r"$\log_{10}\left(|\dot M_{\rm donor}|\right)$")
plt.xlabel("time, $t-t_1$")

plt.subplots_adjust(hspace=0.15)
plt.savefig("paper_figures/separation_mdot_time_log.pdf",bbox_inches='tight')
plt.close()



t1 = ou.get_t1(orb)
tmin = -20
tmax = 2

plt.figure(figsize=(6,6.))

plt.subplot(211)
plt.plot(orb['time']-t1,orb['sep'],'C0-',lw=2,label='simulation')
plt.ylabel(r'separation')
plt.xlim(tmin,tmax)
plt.ylim(0.5,2.06)
plt.grid()
plt.xticks(visible=False)

plt.subplot(212)
plt.plot(md['time']-t1,md['mdot'],'C0-',lw=2 ,label='')
#plt.yscale('log')
plt.xlim(tmin,tmax)
#plt.ylim(1.e-3,1.e-1)
plt.grid()
plt.ylabel(r"$|\dot M_{\rm donor}|$")
plt.xlabel("time, $t-t_1$")
#plt.xticks(visible=False)

plt.subplots_adjust(hspace=0.15)
plt.savefig("paper_figures/separation_mdot_time.pdf",bbox_inches='tight')
plt.close()
