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


mycm = LinearSegmentedColormap.from_list("mycm", ["Black","DodgerBlue"])

thslice=np.pi/2.
phslice=0.0

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
    plt.plot(d['x1v'], np.log10( d["rho"][0,0,:] ),
            color=mycm(float(i)/(nfiles-1)),
            lw=1.3,label=labels[i])
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log \ \rho$ [$M_1/R_1^3$]')
    plt.legend(loc='lower left',fontsize=12,frameon=True)
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(323)
    plt.plot(d['x1v'], np.log10( d["press"][0,0,:] ) ,
            color=mycm(float(i)/(nfiles-1)),
            lw=1.3,label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log \ P$ [$GM_1^2/R_1^4$]')
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(325)
    plt.plot(d['x1v'], d["vel3"][0,0,:],
            color=mycm(float(i)/(nfiles-1)),
            lw=2,
            label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$v_{\phi}$ [$(GM_1/R_1)^{1/2}$]')
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
    plt.plot(d['x1v'], np.log10( d["rho"][0,0,:] ),
            color=mycm(float(i)/(nfiles-1)),
            lw=1.3,label=labels[i])
    plt.xlim(rmin,rmax)
    #plt.ylabel(r'$\log \ \rho$')
    plt.legend(loc='lower left',fontsize=12,frameon=True)
    plt.grid()
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    plt.subplot(324)
    plt.plot(d['x1v'], np.log10( d["press"][0,0,:] ) ,
            color=mycm(float(i)/(nfiles-1)),
            lw=1.3,label='')
    plt.xlim(rmin,rmax)
    #plt.ylabel(r'$\log \ P$')
    plt.grid()
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    plt.subplot(326)
    plt.plot(d['x1v'], d["vel3"][0,0,:],
            color=mycm(float(i)/(nfiles-1)),
            lw=2,
            label='')
    plt.xlim(rmin,rmax)
    #plt.ylabel(r'$v_{\phi}$')
    plt.grid()
    plt.xlabel(r'$r/R_1$')
    plt.yticks(visible=False)
    
plt.subplots_adjust(wspace=0.15,hspace=0.0)
#plt.show()
plt.savefig(plot_name,bbox_inches='tight')
plt.close()




##### 1D & 3D COMPARISON ##########

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
    plt.ylabel(r'$\log \ \rho$')
    plt.legend(loc='lower left',fontsize=12,frameon=True)
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(312)
    plt.plot(d['x1v'], np.log10( d["press"][0,0,:] ) ,
            color=color,ls=linestyle,
            lw=1.3,label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log \ P$')
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
    plt.ylabel(r'$\log \ \rho$')
    plt.legend(loc='lower left',fontsize=12,frameon=True)
    plt.grid()
    plt.xticks(visible=False)
    
    plt.subplot(312)
    plt.plot(d['x1v'], np.log10( np.sum(d["press"][:,0,:],axis=0)/Nphi ) ,
            color=color,ls=linestyle,
            lw=1.3,label='')
    plt.xlim(rmin,rmax)
    plt.ylabel(r'$\log \ P$')
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
