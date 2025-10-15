import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
import OrbitAnalysisUtils as ou
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.colors as colors
import argparse
import seaborn as sns
import pyshtools

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



###### SIMULATION PARAMS   #########
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

G=1
file_list = sorted(glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf"))
file_list = file_list[30::]
print file_list

mylevel=2

####################################

def get_midplane_theta(myfile,level=0):
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)

    # get closest to midplane value
    return dblank['x2v'][ np.argmin(np.abs(dblank['x2v']-np.pi/2.) ) ]

def shift_phi_val(phiarr,valarr,theta_rot):
    """ 2D arrays of phi, value are shifted by theta_rot and re-spliced, returned as 2D """
    print theta_rot
    phases= phiarr + theta_rot
    val =  valarr
    if theta_rot<0:
        #print "theta_rot < 0"
        phases_shift = phases[phases<-np.pi].copy()
        phases_noshift = phases[phases>=-np.pi].copy()
        phases_shift += 2*np.pi
        phases_final = np.concatenate([phases_noshift,phases_shift]).reshape(phiarr.shape)

        val_shift = val[phases<-np.pi].copy()
        val_noshift = val[phases>=-np.pi].copy()
        val_final =  np.concatenate([val_noshift,val_shift]).reshape(phiarr.shape)



    else:
        #print "theta_rot >= 0"
        phases_shift = phases[phases>np.pi].copy()
        phases_noshift = phases[phases<=np.pi].copy()
        phases_shift -= 2*np.pi
        phases_final = np.concatenate([phases_shift,phases_noshift]).reshape(phiarr.shape)

        val_shift = val[phases>np.pi].copy()
        val_noshift = val[phases<=np.pi].copy()
        val_final =  np.concatenate([val_shift,val_noshift]).reshape(phiarr.shape)


    return phases_final,val_final






orb = ou.read_trackfile(base_dir+"pm_trackfile.dat",m1=m1,m2=m2)

print "ORB: ... ", orb.colnames

hst = ascii.read(base_dir+"HSE.hst",
                names=['time','dt','mass','1-mom','2-mom','3-mom','1-KE','2-KE','3-KE','tot-E','mxOmegaEnv','mEnv'])
print "\nHSE: ...", hst.colnames

mg = hst['mass'][0]
print mg+m1

t1=ou.get_t1(orb)
print t1



spec_time=[]
vrms_time=[]


for i,myfile in enumerate(file_list):
    dblank=ar.athdf(myfile,level=mylevel,quantities=[],subsample=True)
    t=dblank['Time']
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
    sma = np.interp(t,orb['time'],orb['sep'])
    Omega = np.interp(t,orb['time'],orb['vmag'])/sma 
    print sma,Omega
    theta_rot = -np.arctan2(y2,x2)
    
    # OVERPLOT THE ROCHE LOBE
    #b = ou.makebinary(1.0,m2,sma)
    #phi = b.get_phi_function()
    #xL,phiL = b.get_xL_phiL()
    #myr= 1.0 #xL[2]-b.x1 #np.abs( xL[0]-b.x1 ) #xL[2]-b.x1
    
    myr = 1.0
    
    
    
    print "r=",myr
    
    x1slicevalue = dblank['x1v'][np.argmin( np.abs(dblank['x1v']-myr) )]
    print "Slicing at x1=",x1slicevalue
    
    x2slicevalue = get_midplane_theta(myfile,level=mylevel)


    # read the midplane data
    
    d=ou.read_data(myfile,orb,m1,m2,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,get_energy=False,
                   x2_min=x2slicevalue,x2_max=x2slicevalue,
                   profile_file=base_dir+"hse_profile.dat")
    
    print "t=",d["Time"]
    
    d['vmag_com'] = np.sqrt((d['vx']-vcom[0])**2 + (d['vy']-vcom[1])**2 + (d['vz']-vcom[2])**2 )
    # ROTATE POSITIONS, VELOCITIES TO COROTATING FRAME
    xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
    yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)
    
    #vxrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.cos(theta_rot) - (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.sin(theta_rot)
    #vyrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.sin(theta_rot) + (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.cos(theta_rot)
    #phicom = np.arctan2(yrot,xrot)
    #vphi = - Omega*np.sqrt(xrot**2 + yrot**2)
    #vxrotC = vxrot - vphi*np.sin(phicom)
    #vyrotC = vyrot + vphi*np.cos(phicom)
    x2rot = (x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot)
    y2rot = (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot)
    


    dr =ou.read_data(myfile,orb,m1,m2,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,get_energy=False,
                     x1_min=x1slicevalue,x1_max=x1slicevalue,
                     profile_file=base_dir+"hse_profile.dat")
    

    #### NOW MAKE THE FIGURE #####
    plt.figure(figsize=(8,8))
    
    plt.subplot2grid((3, 3), (0, 0), colspan=3,rowspan=2)
    im=plt.pcolormesh(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane(np.log10(d['rho'][:,0,:] ) ),
        cmap=plt.cm.magma,rasterized=True,
        vmin=-8,vmax=0)
    
    plt.plot(x2rot,y2rot,'w*',markersize=3)
    
    
    plt.annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-2,1.9),color='w',fontsize='large')
    plt.annotate(r"$a=$"+str(np.round(sma,decimals=2)),(-2,1.65),color='w',fontsize='large')
    plt.axis('equal')
    plt.xlim(-2,3)
    plt.ylim(-1.5,1.5)
    cb = plt.colorbar(im,label=r"$\log_{10} \left(\rho \right)$")
    #plt.show()
    plt.xlabel("$x'$ (donor radii)")
    plt.ylabel("$y'$ (donor radii)")
    plt.tight_layout()
    
    
    
    plt.subplot2grid((3, 3), (2, 0), projection="mollweide",colspan=2)
    
    phiarr,v1arr=shift_phi_val(dr['gx3v'][:,:,0],dr['vel1'][:,:,0],theta_rot)
    #plt.subplot(111, projection="mollweide")
    plt.pcolormesh(phiarr,
                   dr['gx2v'][:,:,0]-np.pi/2,
                   v1arr,
                   cmap=plt.cm.RdBu,
                   vmin=-0.25,vmax=0.25,rasterized=True)
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    cb=plt.colorbar(label=r"radial velocity $(r=1)$",orientation='vertical')
    cb.solids.set_rasterized(True) 
    plt.grid(color='C4')
    
    
    ax3 = plt.subplot2grid((3, 3), (2, 2))
    gd=pyshtools.SHGrid.from_array(dr['vel1'][::2,::1,0].T)
    clm=gd.expand(lmax_calc=10)
    plt.plot(clm.degrees(),clm.spectrum(),lw=2)
    plt.axvline(clm.degrees()[np.argmax(clm.spectrum())],ls='--',color='GoldenRod',lw=2)
    plt.xticks([2,4,6,8,10])
    plt.xlim(2,10)
    plt.ylim(1.e-5,0.1)
    plt.yscale('log')
    plt.grid()
    plt.ylabel(r'power, $S_{ff}(l)$')
    plt.xlabel('degree, $l$')
    plt.tight_layout()

    plt.savefig(output_dir+"modes_density_velocity_spectrum_"+str(i)+".png",dpi=150)
    plt.savefig(output_dir+"modes_density_velocity_spectrum_"+str(i)+".pdf",dpi=150)
    plt.close()


    ## ADD THE SPECTRAL DATA TO THE ARRAY
    entry=[t,sma,Omega]
    for i in range(len(clm.degrees())):
        entry.append(clm.spectrum()[i] )
    spec_time.append(entry)


    ## ADD THE RMS VELOCITIES
    select = ((dr['gx2v']>np.pi/2. - np.pi/32.)&(dr['gx2v']<np.pi/2. + np.pi/32.))
    
    N = len(dr['vel1'][select])
    rms = np.sqrt( np.sum(dr['vel1'][select]**2 ) / N )
    
    rms_m = np.sqrt( np.sum(dr['rho'][select]*dr['dvol'][select] * dr['vel1'][select]**2 ) / 
                     np.sum(dr['rho'][select]*dr['dvol'][select]))

    vrms_time.append([t,sma,rms,rms_m])
    


    
namelist = ['t','sep','Omega']
for i in range(len(clm.degrees())):
    namelist.append( "l"+clm.degrees()[i].astype(str) )
spec_time_table = Table(np.array(spec_time),names=namelist)
ascii.write(spec_time_table,output=output_dir+"mode_spec_time.dat")

vrms_time_table = Table(np.array(vrms_time),names=['t','sep','Vrms','Vrms_m'])
ascii.write(vrms_time_table,output=output_dir+"mode_vrms_eq.dat")
