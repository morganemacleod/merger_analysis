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


def orbit_convergence_plot(base_dir,files,labels,
                           use_t1=False,mycm=plt.cm.viridis):
    
    plt.figure(figsize=(11,5))
    for i,fn in enumerate(files):
        orb = ou.read_trackfile(99,99,base_dir+fn)  # m1,m2 not needed for our analysis
        t1 = ou.get_t1(orb)
    
        plt.subplot(121)
        if use_t1:
            plt.plot(orb['time']-t1,orb['sep'],label=labels[i],
                     lw=2,color=mycm(float(i)/(len(files)-1)))
            plt.xlim(-50,3)
            plt.xlabel('$t-t_1 \ \ [(R_1^3/ G M_1)^{1/2}]$')
        else:
            plt.plot(orb['time'],orb['sep'],label=labels[i],
                    lw=2,color=mycm(float(i)/(len(files)-1)))
            plt.xlabel('$t \ \ [(R_1^3/ G M_1)^{1/2}]$')
            
        plt.grid()
        plt.ylabel('separation $[R_1]$')
        plt.legend(loc='lower left',frameon=True)
    
        plt.subplot(122)
        plt.plot(orb['sep'],orb['ltz'],
                lw=2,color=mycm(float(i)/(len(files)-1)))
        plt.grid()
        plt.xlabel('separation $[R_1]$')
        plt.ylabel(r'$L_{\rm tot,z} \ \ [M_1(GM_1R_1)^{1/2}]$')
    
    plt.subplots_adjust(wspace=0.3)


#### SEPARATION #####
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a185_rin03_r01.dat",
         "pm_trackfile_a19_rin03_r01.dat",
         "pm_trackfile_a2_rin03_r01.dat",
        "pm_trackfile_a206_rin03_r01.dat"]

labels = ['$a_0=1.85R_1$','$a_0=1.9R_1$','$a_0=2.0R_1$','$a_0=2.06R_1$']

mycm = LinearSegmentedColormap.from_list("mycm", ["Black","DodgerBlue"])

orbit_convergence_plot(base_dir,files,labels,use_t1=True, mycm=mycm)
plt.savefig('paper_figures/orbit_convergence_a0.pdf',bbox_inches='tight')
plt.close()


###### SOFTENING #####
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a19_rin03_r005.dat",
         "pm_trackfile_a19_rin03_r0075.dat",
         "pm_trackfile_a19_rin03_r01.dat",
         "pm_trackfile_a19_rin03_r015.dat",
         "pm_trackfile_a19_rin03_r02.dat"]

labels = [r"$r_{\rm soft,2}=0.05 R_1$",
          r"$r_{\rm soft,2}=0.075 R_1$",
         r"$r_{\rm soft,2}=0.1 R_1$",
         r"$r_{\rm soft,2}=0.15 R_1$",
         r"$r_{\rm soft,2}=0.2 R_1$"]

orbit_convergence_plot(base_dir,files,labels,use_t1=False, mycm=mycm)
plt.savefig('paper_figures/orbit_convergence_softening.pdf',bbox_inches='tight')
plt.close()


###### INNER BC ######

base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a19_rin015_r01.dat",
         "pm_trackfile_a19_rin02_r01.dat",
         "pm_trackfile_a19_rin03_r01.dat",
         "pm_trackfile_a19_rin045_r01.dat",
         "pm_trackfile_a19_rin06_r01.dat"]

labels = [r"$r_{\rm in}=0.15 R_1$",
          r"$r_{\rm in}=0.20 R_1$",
         r"$r_{\rm in}=0.30 R_1$",
         r"$r_{\rm in}=0.45 R_1$",
         r"$r_{\rm in}=0.60 R_1$"]

orbit_convergence_plot(base_dir,files,labels,use_t1=False, mycm=mycm)
plt.savefig('paper_figures/orbit_convergence_rin.pdf',bbox_inches='tight')
plt.close()


###### RESOLUTION ########
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a19_rin03_r01_n12.dat",
         "pm_trackfile_a19_rin03_r01.dat",
         "pm_trackfile_a19_rin03_r01_n24.dat",
         "pm_trackfile_a19_rin03_r01_n32.dat",]

labels = ["$12^3$",'$16^3$','$24^3$','$32^3$']

orbit_convergence_plot(base_dir,files,labels,use_t1=False, mycm=mycm)
plt.savefig('paper_figures/orbit_convergence_res.pdf',bbox_inches='tight')





###### FCOROT #######
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a206_rin03_r005_n24_fc0.dat",
         "pm_trackfile_a206_rin03_r005_n24_fc05.dat",
         "pm_trackfile_a206_rin03_r005_n24_fc10.dat"
         ]

labels = [r"$\Omega_{\rm env}=0$",
          r'$\Omega_{\rm env}=0.5 \Omega_{\rm orb}$',
          r'$\Omega_{\rm env}= \Omega_{\rm orb}$']

plt.figure()
plt.axhline(2.06,ls='--',lw=1,color='grey')
for i,fn in enumerate(files):
    orb = ou.read_trackfile(99,99,base_dir+fn)  # m1,m2 not needed for our analysis
    t1 = ou.get_t1(orb)
    plt.plot(orb['time']-t1,orb['sep'],label=labels[i],
                     lw=2,color=mycm(float(i)/(len(files)-1)))
    
plt.grid()
plt.ylabel('separation $[R_1]$')
plt.xlabel('$t-t_1 \ \ [(R_1^3/ G M_1)^{1/2}]$')
plt.legend(loc='lower left',frameon=True)
plt.ylim(1.5,)
plt.xlim(-250,)
plt.savefig('paper_figures/orbit_fcorot.pdf',bbox_inches='tight')
plt.close()


plt.figure()

for i,fn in enumerate(files):
    orb = ou.read_trackfile(99,99,base_dir+fn)  # m1,m2 not needed for our analysis
    # number of orbits for decay of angular momentum
    Omega = orb['vmag']/orb['sep']
    dadt = np.gradient(orb['sep'])/np.gradient(orb['time'])
    
    plot_arr = orb['sep']/np.abs(dadt) * Omega/(2.*np.pi)

    plt.plot(orb['sep'],plot_arr,
         label=labels[i],color=mycm(float(i)/(len(files)-1)) )

    plt.legend(loc=0)

plt.yscale('log')
plt.xlabel(r'separation $[R_1]$')
plt.ylabel(r'$N_{\rm decay}$')
plt.ylim(0.3,1e2)
plt.xlim(0.6,2.06)
plt.savefig('paper_figures/ndecay_fcorot.pdf',bbox_inches='tight')
plt.close()
