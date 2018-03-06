# normal stuff
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap


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


##### RES COMPARISON 1D ######
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/1d/"
folders = ['r1536/','r768/','r384/','r192/','r96/']
filename = "HSE.out1.00010.athdf"

plot_name="paper_figures/hse_1D_res_study_midplane.pdf"

dr_o_r96 = 0.05 

# labels 
labels = [r'$dr/r='+str(np.round(dr_o_r96/16,decimals=4))+'$',
          r'$dr/r='+str(np.round(dr_o_r96/8,decimals=4))+'$',
          r'$dr/r='+str(np.round(dr_o_r96/4,decimals=4))+'$',
          r'$dr/r='+str(np.round(dr_o_r96/2,decimals=4))+'$',
          r'$dr/r='+str(np.round(dr_o_r96/1,decimals=4))+'$']

files = [base_dir+folders[0]+filename,
        base_dir+folders[1]+filename,
        base_dir+folders[2]+filename,
        base_dir+folders[3]+filename,
        base_dir+folders[4]+filename]

# get the slice vals
d=ar.athdf(files[0],quantities=[])
thval = np.pi/2
phval = 0.0

thslice = d['x2v'][ np.argmin( np.abs(d['x2v']-thval) ) ]
phslice = d['x3v'][ np.argmin( np.abs(d['x3v']-phval) ) ]
x3ind = 0
x2ind = 0
Omega = 0

mycm = LinearSegmentedColormap.from_list("mycm", ["Black","DodgerBlue"])


rmin = 0.25
rmax = 1.5

plt.figure(figsize=(11,8))
nfiles = len(files)

for i,fn in enumerate(files):
    d=ar.athdf(fn,level=0,
           x1_max=2,
           x2_min=thslice,x2_max=thslice,
           x3_min=phslice,x3_max=phslice)
    
    
    plt.subplot(321)
    plt.plot(d['x1v'], np.log10( d["rho"][x3ind,x2ind,:] ),
            color=mycm(float(i)/(nfiles-1)),
            lw=1.3,label=labels[i])
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log_{10} \left( \rho \right)$')
    plt.legend(loc='lower left',fontsize=12,frameon=True)
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(323)
    plt.plot(d['x1v'], np.log10( d["press"][x3ind,x2ind,:] ) ,
            color=mycm(float(i)/(nfiles-1)),
            lw=1.3,label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log_{10} \left( P \right)$')
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(325)
    plt.plot(d['x1v'], d["vel3"][x3ind,x2ind,:] + Omega*d['x1v']*np.sin(d['x2v']),
            color=mycm(float(i)/(nfiles-1)),
            lw=2,
            label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$v_{\phi}$')
    plt.grid()
    plt.xlabel(r'$r/R_1$')
    
    
# NOW THE ZOOM
rmin = 0.9
rmax = 1.1

nfiles = len(files)

for i,fn in enumerate(files):
    d=ar.athdf(fn,level=0,
           x1_max=2,
           x2_min=thslice,x2_max=thslice,
           x3_min=phslice,x3_max=phslice)
    
    
    plt.subplot(322)
    plt.plot(d['x1v'], np.log10( d["rho"][x3ind,x2ind,:] ),
            color=mycm(float(i)/(nfiles-1)),
            lw=1.3,label=labels[i])
    plt.xlim(rmin,rmax)
    #plt.ylabel(r'$\log \ \rho$')
    plt.legend(loc='lower left',fontsize=12,frameon=True)
    plt.grid()
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    plt.subplot(324)
    plt.plot(d['x1v'], np.log10( d["press"][x3ind,x2ind,:] ) ,
            color=mycm(float(i)/(nfiles-1)),
            lw=1.3,label='')
    plt.xlim(rmin,rmax)
    #plt.ylabel(r'$\log \ P$')
    plt.grid()
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    plt.subplot(326)
    plt.plot(d['x1v'], d["vel3"][x3ind,x2ind,:] + Omega*d['x1v']*np.sin(d['x2v']),
            color=mycm(float(i)/(nfiles-1)),
            lw=2,
            label='')
    plt.xlim(rmin,rmax)
    #plt.ylabel(r'$v_{\phi}$')
    plt.grid()
    plt.xlabel(r'radius')
    plt.yticks(visible=False)
    
plt.subplots_adjust(wspace=0.15,hspace=0.0)
#plt.show()
plt.savefig(plot_name,bbox_inches='tight')
plt.close()




### 1D & 3D

def read_and_plot(fn,thslice,phslice,label,
                 color='k',linestyle='-',
                 level=2):
    d=ar.athdf(fn,level=level,
           x1_max=2,
           x2_min=thslice,x2_max=thslice,
           x3_min=phslice,x3_max=phslice)
    
    
    plt.subplot(311)
    plt.plot(d['x1v'], np.log10( d["rho"][0,0,:] ),
            color=color,ls=linestyle,
            lw=1.3,label=label)
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log_{10} \left( \rho \right)$')
    plt.legend(loc='lower left',fontsize=12,frameon=True)
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(312)
    plt.plot(d['x1v'], np.log10( d["press"][0,0,:] ) ,
            color=color,ls=linestyle,
            lw=1.3,label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log_{10} \left( P \right)$')
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(313)
    plt.plot(d['x1v'], d["vel3"][0,0,:],
            color=color,ls=linestyle,
            lw=2,
            label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$v_{\phi}$')
    plt.grid()
    plt.xlabel(r'$r$')
    
# GET THE SLICE OF CHOICE
def get_thslice_phislice(fn,thval,phval,level=2):
    d=ar.athdf(fn,quantities=[],level=level)

    thslice = d['x2v'][ np.argmin( np.abs(d['x2v']-thval) ) ]
    phslice = d['x3v'][ np.argmin( np.abs(d['x3v']-phval) ) ]

    return thslice,phslice

def read_and_plot_phiavg(fn,thslice,phslice,label,
                 color='k',linestyle='-',
                 level=2):
    d=ar.athdf(fn,level=level,
           x1_max=2,
           x2_min=thslice,x2_max=thslice)
    
    Nphi = len(d["rho"][:,0,0])

    
    
    plt.subplot(311)
    plt.plot(d['x1v'], np.log10( np.sum(d["rho"][:,0,:],axis=0)/Nphi ),
            color=color,ls=linestyle,
            lw=1.3,label=label)
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log_{10} \left( \rho \right)$')
    plt.legend(loc='lower left',fontsize=12,frameon=True)
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(312)
    plt.plot(d['x1v'], np.log10( np.sum(d["press"][:,0,:],axis=0)/Nphi ) ,
            color=color,ls=linestyle,
            lw=1.3,label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log_{10} \left( P \right)$')
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(313)
    plt.plot(d['x1v'], np.sum(d["vel3"][:,0,:],axis=0)/Nphi,
            color=color,ls=linestyle,
            lw=2,
            label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$v_{\phi}$')
    plt.grid()
    plt.xlabel(r'$r$')


### MAKE THE PLOT
threeDfile= glob("/Volumes/DATAVolume/athenaruns/pm_envelope/pole/hse/res24_rot/HSE.out1.000[0-9][0-9].athdf")[-1]
oneDfile = glob("/Volumes/DATAVolume/athenaruns/pm_envelope/1d/r1536/HSE.out1.000[0-9][0-9].athdf")[-1]

rmin = 0.25
rmax = 1.5

plt.figure(figsize=(6,8))

thslice,phslice = get_thslice_phislice(threeDfile,0,0)
read_and_plot(threeDfile,thslice,phslice,r"3D, $\theta="+str(np.round(thslice,decimals=3))+"$" ,level=0)

thslice,phslice = get_thslice_phislice(threeDfile,np.pi/2,0)
read_and_plot(threeDfile,thslice,phslice, r"3D, $\theta="+str(np.round(thslice,decimals=3))+"$",
             color='r',level=2)

thslice,phslice = get_thslice_phislice(oneDfile,np.pi/2,0)
read_and_plot(oneDfile,thslice,phslice, r"1D, $\theta="+str(np.round(thslice,decimals=3))+"$",
             color='b',level=0,linestyle='--')



plt.subplots_adjust(hspace=0)
plt.savefig("paper_figures/hse_1d_3d_comparison.pdf",bbox_inches='tight')
plt.close()



##### ORBIT CONVERGENCE PLOTS

import OrbitAnalysisUtils as ou

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
            plt.xlabel('time, $t-t_1$')
        else:
            plt.plot(orb['time'],orb['sep'],label=labels[i],
                    lw=2,color=mycm(float(i)/(len(files)-1)))
            plt.xlabel('time, $t$')
            
        plt.grid()
        plt.ylabel('separation')
        plt.legend(loc='lower left',frameon=True)
    
        plt.subplot(122)
        plt.plot(orb['sep'],orb['ltz'],
                lw=2,color=mycm(float(i)/(len(files)-1)))
        plt.grid()
        plt.xlabel('separation')
        plt.ylabel(r'total $\hat z$ angular momentum')
    
    plt.subplots_adjust(wspace=0.3)

### Separation
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a185_rin03_r01.dat",
         "pm_trackfile_a19_rin03_r01.dat",
         "pm_trackfile_a2_rin03_r01.dat",
        "pm_trackfile_a206_rin03_r01.dat"]

labels = ['$a_0=1.85$','$a_0=1.9$','$a_0=2.0$','$a_0=2.06$']

mycm = LinearSegmentedColormap.from_list("mycm", ["Black","DodgerBlue"])

orbit_convergence_plot(base_dir,files,labels,use_t1=True, mycm=mycm)
plt.savefig('paper_figures/orbit_convergence_a0.pdf',bbox_inches='tight')

### Softening
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a19_rin03_r005.dat",
         "pm_trackfile_a19_rin03_r0075.dat",
         "pm_trackfile_a19_rin03_r01.dat",
         "pm_trackfile_a19_rin03_r015.dat",
         "pm_trackfile_a19_rin03_r02.dat"]

labels = [r"$r_{\rm soft,2}=0.05$",
          r"$r_{\rm soft,2}=0.075$",
         r"$r_{\rm soft,2}=0.1$",
         r"$r_{\rm soft,2}=0.15$",
         r"$r_{\rm soft,2}=0.2$"]

orbit_convergence_plot(base_dir,files,labels,use_t1=False, mycm=mycm)
plt.savefig('paper_figures/orbit_convergence_softening.pdf',bbox_inches='tight')

### rin
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a19_rin015_r01.dat",
         "pm_trackfile_a19_rin02_r01.dat",
         "pm_trackfile_a19_rin03_r01.dat",
         "pm_trackfile_a19_rin045_r01.dat",
         "pm_trackfile_a19_rin06_r01.dat"]

labels = [r"$r_{\rm in}=0.15$",
          r"$r_{\rm in}=0.20$",
         r"$r_{\rm in}=0.30$",
         r"$r_{\rm in}=0.45$",
         r"$r_{\rm in}=0.60$"]

orbit_convergence_plot(base_dir,files,labels,use_t1=False, mycm=mycm)
plt.savefig('paper_figures/orbit_convergence_rin.pdf',bbox_inches='tight')

### res
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

files = ["pm_trackfile_a19_rin03_r01_n12.dat",
         "pm_trackfile_a19_rin03_r01.dat",
         "pm_trackfile_a19_rin03_r01_n24.dat",
         "pm_trackfile_a19_rin03_r01_n32.dat",]

labels = ["$12^3$",'$16^3$','$24^3$','$32^3$']

orbit_convergence_plot(base_dir,files,labels,use_t1=False, mycm=mycm)
plt.savefig('paper_figures/orbit_convergence_res.pdf',bbox_inches='tight')



#### COROTATION FRACTION
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/trackfiles/"

mycm = LinearSegmentedColormap.from_list("mycm", ["Black","DodgerBlue"])

files = ["pm_trackfile_a206_rin03_r005_n24_fc0.dat",
         "pm_trackfile_a206_rin03_r005_n24_fc05.dat",
         "pm_trackfile_a206_rin03_r005_n24_fc10.dat"
         ]

labels = [r"$\Omega_{\rm env}=0$",
          r'$\Omega_{\rm env}=0.5 \Omega_{\rm orb,0}$',
          r'$\Omega_{\rm env}= \Omega_{\rm orb,0}$']

plt.figure()
plt.axhline(2.06,ls='-',lw=1,color='SkyBlue',zorder=0)
plt.annotate(r"$a_{\rm RL}$",(3,2.03))
for i,fn in enumerate(files):
    orb = ou.read_trackfile(99,99,base_dir+fn)  # m1,m2 not needed for our analysis
    t1 = ou.get_t1(orb)
    plt.plot(orb['time']-t1,orb['sep'],label=labels[i],
                     lw=2,color=mycm(float(i)/(len(files)-1)))
    
plt.grid()
plt.ylabel('separation')
plt.xlabel('time, $t-t_1$')
plt.legend(loc='lower left',frameon=True)
plt.ylim(1.5,)
plt.xlim(-250,)

plt.savefig('paper_figures/fcorot_sep_time.pdf',bbox_inches='tight')

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

plt.axvline(1,color='grey',ls='--',zorder=0)
plt.annotate(r"$R_1$",(1.03,5.5))

plt.yscale('log')
plt.xlabel(r'separation')
plt.ylabel(r'$N_{\rm decay}$')
plt.ylim(0.3,1e2)
plt.xlim(0.6,2.06)
plt.savefig('paper_figures/fcorot_ndecay.pdf',bbox_inches='tight')
