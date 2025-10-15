import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import athena_read as ar
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
import OrbitAnalysisUtils as ou
from scipy.misc import derivative
from scipy.integrate import odeint
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import ImageGrid



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



# derivatives function 
def derivs(vec,time):
    x = vec[0]
    y = vec[1]
    z = vec[2]
    r = np.array([x,y,z])
    vx = vec[3]
    vy = vec[4]
    vz = vec[5]
    v = np.array([vx,vy,vz])
    
    dist1 = np.sqrt( (x-b.x1)**2 + y**2 + z**2)
    dist2 = np.sqrt( (x-b.x2)**2 + y**2 + z**2)
    
    centrifugal = - np.cross(b.omega_vec,np.cross(b.omega_vec,r) )
    coriolis = -2*np.cross(b.omega_vec,v)
    
    dxdt = vx
    dydt = vy
    dzdt = vz
    dvxdt = - m1/dist1**3*(x-b.x1) - m2*ou.fspline(dist2,0.1)*(x-b.x2) + centrifugal[0] + coriolis[0]
    dvydt = - m1/dist1**3*y        - m2*ou.fspline(dist2,0.1)*y        + centrifugal[1] + coriolis[1]
    dvzdt = - m1/dist1**3*z        - m2*ou.fspline(dist2,0.1)*z        + centrifugal[2] + coriolis[2]
    
    derivsarray = np.array([dxdt,dydt,dzdt,dvxdt,dvydt,dvzdt])
    
    return derivsarray


def binary_balistic(initialvalues,tlim=1,ntimes=1001):
    """ returns an ascii table of integrated trajectory given input values"""
    times = np.linspace(0,tlim,ntimes)
    ode_result = odeint(derivs, initialvalues, times)
    # Convert the result to a table
    column_names = ['x','y','z','vx','vy','vz']
    sol = Table(ode_result,names=column_names)
    # add the times column
    sol['t'] = times
    return sol




def interpolate_vxC_vyC(interpoints):
    points = np.array([xrot.T,yrot.T]).T
    points = points.reshape(len(points[0])*len(points[:,0]),2)
    vxd = vxrotC.reshape(len(vxrotC[0])*len(vxrotC[:,0]))
    vyd = vyrotC.reshape(len(vyrotC[0])*len(vyrotC[:,0]))
    
    vxi = griddata(points,vxd, interpoints,method='nearest')
    vyi = griddata(points,vyd, interpoints,method='nearest')

    return np.array([vxi.T,vyi.T]).T



###### SIMULATION PARAMS   #########
base_dir = "/Volumes/DATAVolume/athenaruns/pm_envelope/pole/syncRL/fcorot_series/a206_res24_fc10/" 

output_dir = "paper_figures/"

m1 = 0.410103
m2 = 0.3
G=1

full_file_list = glob(base_dir+"HSE.out1.00[0-9][0-9][0-9].athdf")

file_list = [full_file_list[-322],
             full_file_list[-13]]




mylevel=2


####################################

def get_midplane_theta(myfile,level=0):
    dblank=ar.athdf(myfile,level=level,quantities=[],subsample=True)

    # get closest to midplane value
    return dblank['x2v'][ np.argmin(np.abs(dblank['x2v']-np.pi/2.) ) ]



orb = ou.read_trackfile(m1,m2,base_dir+"pm_trackfile.dat")
t1=ou.get_t1(orb)
print "t1=",t1





x2slicevalue=get_midplane_theta(file_list[0],level=mylevel)
print "Slicing at x2=",x2slicevalue


fig = plt.figure(1,figsize=(10,9))
nrows = 2
ncols = 3
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(nrows,ncols),  # creates 2x2 grid of axes
                     axes_pad=0.33,  # pad between axes in inch.
                     cbar_mode='edge',
                     cbar_location='top',
                     cbar_size="2%",
                     cbar_pad="2%")
    
for i,myfile in enumerate(file_list):
    print i, myfile 
        
    # read the data
    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x2_min=x2slicevalue,x2_max=x2slicevalue)
    t = d['Time']
    print "t=",t
        
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
        
    theta_rot = -np.arctan2(y2,x2)
    
    sma = np.interp(t,orb['time'],orb['sep'])
    Omega = np.interp(t,orb['time'],orb['vmag'])/sma 
    print sma,Omega

    # ROTATE POSITIONS, VELOCITIES TO COROTATING FRAME
    xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
    yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)

    vxrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.cos(theta_rot) - (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.sin(theta_rot)
    vyrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.sin(theta_rot) + (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.cos(theta_rot)


    phicom = np.arctan2(yrot,xrot)

    vphi = - Omega*np.sqrt(xrot**2 + yrot**2)
    vxrotC = vxrot - vphi*np.sin(phicom)
    vyrotC = vyrot + vphi*np.cos(phicom)
    
    ##########
    #  DENSITY
    #########
    ind = 3*i+0
    im=grid[ind].pcolormesh(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane(np.log10(d['rho'][:,len(d['x2v'])/2,:]) ),
        cmap=plt.cm.magma,
        vmin=-8,vmax=0,rasterized=True)
    
    x2rot = (x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot)
    y2rot = (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot)
    x1rot = (-rcom[0])*np.cos(theta_rot)-(-rcom[1])*np.sin(theta_rot)
    y1rot = (-rcom[0])*np.sin(theta_rot)+(-rcom[1])*np.cos(theta_rot)
    grid[ind].plot(x2rot,y2rot,'w*',markersize=3)
    grid[ind].plot(x1rot,y1rot,'b*',markersize=3)
    grid[ind].plot(0,0,'kx',markersize=3)
    
    grid[ind].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(-0.3,1.2),color='w',fontsize='large')
        
        
    grid[ind].set_xlim(-0.5,2.5)
    grid[ind].set_ylim(-1.5,1.5)
    grid[ind].set_xlabel(r"$x'$")
    grid[ind].set_ylabel(r"$y'$")
    cb = grid.cbar_axes[ind].colorbar(im)
    cb.set_label_text(r"$\log_{10} \left( \rho \right)$")
        
        
    # OVERPLOT THE ROCHE LOBE
    b = ou.makebinary(1.0,m2,sma)
    phi = b.get_phi_function()

    # PLOT THE CONTOUR
    xL,phiL = b.get_xL_phiL()
    xl= np.linspace(-0.5,2.5,301)
    yl = np.linspace(-1.5,1.5,301)
    xx,yy = np.meshgrid(xl,yl)
    grid[ind].contour(xx,yy,phi(xx,yy,0),levels=np.sort(phiL),colors='C4',linestyles='-',linewidths=0.5,rasterized=True)


        
        
        
        
        
    ##########
    #  Pressure
    #########
    ind = 3*i+1
    im=grid[ind].pcolormesh(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane(np.log10(d['press'][:,len(d['x2v'])/2,:]) ),
        cmap=plt.cm.viridis,
        vmin=-8,vmax=0,rasterized=True)
    
    grid[ind].plot(x2rot,y2rot,'w*',markersize=3)
    grid[ind].plot(x1rot,y1rot,'b*',markersize=3)
    grid[ind].contour(xx,yy,phi(xx,yy,0),levels=np.sort(phiL),colors='C4',linestyles='-',linewidths=0.5,rasterized=True)
    
    grid[ind].set_xlim(-0.5,2.5)
    grid[ind].set_ylim(-1.5,1.5)
    grid[ind].set_xlabel(r"$x'$")
    grid[ind].set_ylabel(r"$y'$")
    cb = grid.cbar_axes[ind].colorbar(im)
    cb.set_label_text(r"$\log_{10} \left( P \right)$")
        
        
    ##########
    #  Entropy
    #########
    ind = 3*i+2
    im=grid[ind].pcolormesh(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane(np.log(d['press'][:,len(d['x2v'])/2,:]/d['rho'][:,len(d['x2v'])/2,:]**(5./3.)) ),
        cmap=plt.cm.RdYlBu_r,
        vmin=-0.5,vmax=5.5,rasterized=True)
    
    grid[ind].plot(x2rot,y2rot,'w*',markersize=3)
    grid[ind].plot(x1rot,y1rot,'b*',markersize=3)
    
    skip = 16
    grid[ind].quiver(xrot[::skip,::skip],
           yrot[::skip,::skip],
           vxrotC[::skip,::skip],
           vyrotC[::skip,::skip],
           scale=4,scale_units='xy'
           )
        
        
    grid[ind].set_xlim(-0.5,2.5)
    grid[ind].set_ylim(-1.5,1.5)
    grid[ind].set_xlabel(r"$x'$")
    grid[ind].set_ylabel(r"$y'$")
    cb = grid.cbar_axes[ind].colorbar(im)
    cb.set_label_text(r"$\ln \left(P/\rho^\gamma \right)$")
    
plt.savefig(output_dir+"flow_zoom_secondary.pdf",bbox_inches='tight',dpi=150)
plt.clf()
















###### SIMULATION PARAMS   #########
file_list = [full_file_list[-322],
             full_file_list[-62],
             full_file_list[-32],
             full_file_list[-13]]

####################################



fig = plt.figure(1,figsize=(11,10))
nrows = 2
ncols = 2
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(nrows,ncols),  # creates 2x2 grid of axes
                     axes_pad=0.33,  # pad between axes in inch.
                     cbar_mode='edge',
                     cbar_location='right',
                     cbar_size="2%",
                     cbar_pad="2%")
    
for i,myfile in enumerate(file_list):
    print i, myfile 
        
    # read the data
    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x2_min=x2slicevalue,x2_max=x2slicevalue)
    t = d['Time']
        
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
        
    theta_rot = -np.arctan2(y2,x2)
    
    sma = np.interp(t,orb['time'],orb['sep'])
    Omega = np.interp(t,orb['time'],orb['vmag'])/sma 
    print sma,Omega

    # ROTATE POSITIONS, VELOCITIES TO COROTATING FRAME
    xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
    yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)

    vxrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.cos(theta_rot) - (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.sin(theta_rot)
    vyrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.sin(theta_rot) + (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.cos(theta_rot)


    phicom = np.arctan2(yrot,xrot)

    vphi = - Omega*np.sqrt(xrot**2 + yrot**2)
    vxrotC = vxrot - vphi*np.sin(phicom)
    vyrotC = vyrot + vphi*np.cos(phicom)
    
    x2rot = (x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot)
    y2rot = (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot)
    x1rot = (-rcom[0])*np.cos(theta_rot)-(-rcom[1])*np.sin(theta_rot)
    y1rot = (-rcom[0])*np.sin(theta_rot)+(-rcom[1])*np.cos(theta_rot) 



        
    ##########
    #  Entropy
    #########
    ind = i
    im=grid[ind].pcolormesh(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane(np.log(d['press'][:,len(d['x2v'])/2,:]/d['rho'][:,len(d['x2v'])/2,:]**(5./3.)) ),
        cmap=plt.cm.RdYlBu_r,
        vmin=-0.5,vmax=5.5,rasterized=True)
    
    grid[ind].plot(x2rot,y2rot,'w*',markersize=3)
    grid[ind].plot(x1rot,y1rot,'b*',markersize=3)    
    
    
    # OVERPLOT THE ROCHE LOBE
    b = ou.makebinary(1.0,m2,sma)
    phi = b.get_phi_function()

    # PLOT THE CONTOUR
    xL,phiL = b.get_xL_phiL()
    xl= np.linspace(-3,3,301)
    yl = np.linspace(-3,3,301)
    xx,yy = np.meshgrid(xl,yl)
    grid[ind].contour(xx,yy,phi(xx,yy,0),levels=np.sort(phiL),colors='C4',linestyles='-',lw=1,rasterized=True)
    
    
    
    # CONSTRUCT INITIAL VALUES
    width = 0.06
    npoints = 20
    x0 = xL[1] + np.random.normal(0,width/2,npoints) + width
    y0 = 0 + np.random.normal(0,width,npoints)
    z0 = np.zeros_like(y0)

    # interpolate 
    interpoints = np.array([x0.T,y0.T]).T
    vxvy = interpolate_vxC_vyC(interpoints)

    print vxvy

    vx0 = vxvy[:,0]
    vy0 = vxvy[:,1]
    vz0 = np.zeros_like(y0)


    initialvalues = np.array([x0.T,y0.T,z0.T,vx0.T,vy0.T,vz0.T]).T


    
    

    skip = 15
    grid[ind].quiver(xrot[::skip,::skip],
           yrot[::skip,::skip],
           vxrotC[::skip,::skip],
           vyrotC[::skip,::skip],
           scale=5,scale_units='xy'
           )
        
    # PLOT THE INTEGRATED TRAJ
    for j,init in enumerate(initialvalues):
        sol = binary_balistic(init,tlim=10,ntimes=1001)
        # PLOT THE SOLUTION
        grid[ind].plot(sol['x'],sol['y'],'-',lw=1,color='skyblue',rasterized=True)




    grid[ind].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(0.9,1.2),color='k',fontsize='large',backgroundcolor='FloralWhite')
    grid[ind].set_xlim(-0.5,2.5)
    grid[ind].set_ylim(-1.5,1.5)
    grid[ind].set_xlabel(r"$x'$")
    grid[ind].set_ylabel(r"$y'$")
    cb = grid.cbar_axes[ind].colorbar(im)
    cb.set_label_text(r"$\ln \left(P/\rho^\gamma \right)$")
    
plt.savefig(output_dir+"ballistic_zoom_secondary.pdf",bbox_inches='tight',dpi=150)
plt.clf()








#######################
#  OMEGA
#######################

fcorotation=1.0

sma0 = np.interp(0,orb['time'],orb['sep'])
Omega0 = fcorotation*np.interp(0,orb['time'],orb['vmag'])/sma0 

fig = plt.figure(1,figsize=(7,6))
nrows = 2
ncols = 2
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(nrows,ncols),  # creates 2x2 grid of axes
                     axes_pad=0.2,  # pad between axes in inch.
                     cbar_mode='single',
                     cbar_location='right',
                     cbar_size="2%",
                     cbar_pad="2%")
    
for i,myfile in enumerate(file_list):
    print i, myfile 
        
    # read the data
    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x2_min=x2slicevalue,x2_max=x2slicevalue)
    t = d['Time']
        
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
        
    theta_rot = -np.arctan2(y2,x2)
    
    sma = np.interp(t,orb['time'],orb['sep'])
    Omega = np.interp(t,orb['time'],orb['vmag'])/sma 
    print sma,Omega

    # ROTATE POSITIONS, VELOCITIES TO COROTATING FRAME
    xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
    yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)

    vxrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.cos(theta_rot) - (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.sin(theta_rot)
    vyrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.sin(theta_rot) + (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.cos(theta_rot)


    phicom = np.arctan2(yrot,xrot)

    vphi = - Omega*np.sqrt(xrot**2 + yrot**2)
    vxrotC = vxrot - vphi*np.sin(phicom)
    vyrotC = vyrot + vphi*np.cos(phicom)
    
    x2rot = (x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot)
    y2rot = (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot)
    x1rot = (-rcom[0])*np.cos(theta_rot)-(-rcom[1])*np.sin(theta_rot)
    y1rot = (-rcom[0])*np.sin(theta_rot)+(-rcom[1])*np.cos(theta_rot)
 
    OmegaPlot = d['vel3'][:,len(d['x2v'])/2,:]/d['gx1v'][:,len(d['x2v'])/2,:] - Omega


        
    ##########
    # Omega
    #########
    ind = i
    im=grid[ind].pcolormesh(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane( OmegaPlot ),
        cmap=plt.cm.PuOr,
        vmin=-1,vmax=1,rasterized=True)
        

    grid[ind].contour(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane( np.log10(d['rho'][:,len(d['x2v'])/2,:]) ),
        cmap=plt.cm.bone_r,linestyles='-',levels=np.linspace(-8,0,17),rasterized=True)

   
    grid[ind].plot(x2rot,y2rot,'w*',markersize=3)
    grid[ind].plot(x1rot,y1rot,'b*',markersize=3)
    
    grid[ind].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(0.5,1.2),color='k',fontsize='large',backgroundcolor='FloralWhite')
    grid[ind].set_xlim(-1.5,2.5)
    grid[ind].set_ylim(-1.5,1.5)
    grid[ind].set_xlabel(r"$x'$")
    grid[ind].set_ylabel(r"$y'$")
    cb = grid.cbar_axes[ind].colorbar(im)
    cb.set_label_text(r"$\Omega-\Omega_{\rm orb}$")
    
plt.savefig(output_dir+"omega_inst_zoom_secondary.pdf",bbox_inches='tight',dpi=150)
plt.clf()



fig = plt.figure(1,figsize=(7,6))
nrows = 2
ncols = 2
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(nrows,ncols),  # creates 2x2 grid of axes
                     axes_pad=0.2,  # pad between axes in inch.
                     cbar_mode='single',
                     cbar_location='right',
                     cbar_size="2%",
                     cbar_pad="2%")
    
for i,myfile in enumerate(file_list):
    print i, myfile 
        
    # read the data
    d = ou.read_data(myfile,orb,m1,m2,G=1,rsoft2=0.05,level=mylevel,get_cartesian=True,get_torque=False,
                         x1_max=7.5,x2_min=x2slicevalue,x2_max=x2slicevalue)
    t = d['Time']
        
    rcom,vcom = ou.rcom_vcom(orb,t)
    x2,y2,z2 = ou.pos_secondary(orb,t)
        
    theta_rot = -np.arctan2(y2,x2)
    
    sma = np.interp(t,orb['time'],orb['sep'])
    Omega = np.interp(t,orb['time'],orb['vmag'])/sma 
    print sma,Omega

    # ROTATE POSITIONS, VELOCITIES TO COROTATING FRAME
    xrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.cos(theta_rot) - (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.sin(theta_rot)
    yrot = (d['x'][:,len(d['x2v'])/2,:]-rcom[0])*np.sin(theta_rot) + (d['y'][:,len(d['x2v'])/2,:]-rcom[1])*np.cos(theta_rot)

    vxrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.cos(theta_rot) - (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.sin(theta_rot)
    vyrot = (d['vx'][:,len(d['x2v'])/2,:]-vcom[0])*np.sin(theta_rot) + (d['vy'][:,len(d['x2v'])/2,:]-vcom[1])*np.cos(theta_rot)


    phicom = np.arctan2(yrot,xrot)

    vphi = - Omega*np.sqrt(xrot**2 + yrot**2)
    vxrotC = vxrot - vphi*np.sin(phicom)
    vyrotC = vyrot + vphi*np.cos(phicom)
    
    x2rot = (x2-rcom[0])*np.cos(theta_rot)-(y2-rcom[1])*np.sin(theta_rot)
    y2rot = (x2-rcom[0])*np.sin(theta_rot)+(y2-rcom[1])*np.cos(theta_rot)
 
    OmegaPlot = d['vel3'][:,len(d['x2v'])/2,:]/d['gx1v'][:,len(d['x2v'])/2,:] - Omega0


        
    ##########
    #  Omega
    #########
    ind = i
    im=grid[ind].pcolormesh(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane( OmegaPlot ),
        cmap=plt.cm.PuOr,
        vmin=-1,vmax=1,rasterized=True)
    
    grid[ind].plot(x2rot,y2rot,'w*',markersize=3)
        

    grid[ind].contour(
        ou.get_plot_array_midplane(xrot),
        ou.get_plot_array_midplane(yrot),
        ou.get_plot_array_midplane( np.log10(d['rho'][:,len(d['x2v'])/2,:]) ),
        cmap=plt.cm.bone_r,linestyles='-',levels=np.linspace(-8,0,17),rasterized=True)
    
    grid[ind].annotate(r"$t-t_1=$"+str(np.round(t-t1,decimals=2)),(0.5,1.2),color='k',fontsize='large',backgroundcolor='FloralWhite')
    grid[ind].set_xlim(-1.5,2.5)
    grid[ind].set_ylim(-1.5,1.5)
    grid[ind].set_xlabel(r"$x'$")
    grid[ind].set_ylabel(r"$y'$")
    cb = grid.cbar_axes[ind].colorbar(im)
    cb.set_label_text(r"$\Omega-\Omega_{\rm orb,0}$")
    
plt.savefig(output_dir+"omega_zoom_secondary.pdf",bbox_inches='tight',dpi=150)
plt.clf()
