import matplotlib as mpl
mpl.use('agg')
import numpy as np
import matplotlib.pyplot as plt
import athena_read as ar
import OrbitAnalysisUtils as ou

from Constants import Constants
#import seaborn as sns
#from astropy.io import ascii
c=Constants()

from skimage import measure
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d

import deepdish as dd
from astropy.table import Table
from glob import glob
import argparse
import multiprocessing as mp
import os

X = 0.73
Z = 0.02
def kappa(rho,T, returnparts=False):
    """From G. Knapp A403 Princeton Course Notes"""
    Kes = 0.2*(1+X)*rho**0 *T**0
    Ke = 0.2*(1+X)*(1+2.7e11*rho/T**2)**(-1) *(1+ (T/4.5e8)**0.86)**(-1)
    Kk = 4.e25*(1+X)*(Z+1.e-3)*rho/T**3.5
    Khm = 1.1e-25*Z**0.5 *rho**0.5 * T**7.7
    Km = 0.1*Z*rho**0 * T**0
    Krad = Km + (Khm**-1 + (Ke+Kk)**-1 )**-1
    if returnparts:
        return Krad, [Kes,Ke,Kk,Khm,Km]
    else:
        return Krad

def get_interp_function(d,var,method='linear'):
    """
    MM: Use RegularGridInterpolator to pass data to interpolating function for a given variable
    Parameters
    -----------
    d : dict
       athena data dict from read_data
    var: str
       name of variable to be interpolated
       
    Returns
    --------
    var_interp: an interpolating function that can be called with a tuple (phi,theta,r)
    """
    dph = np.gradient(d['x3v'])[0]
    two_pi = ( (d['x3v'][-1]-d['x3v'][0]+dph) /(2*np.pi) > 0.99 ) # boolean to determine if spans 2pi in phi
    
    if two_pi:
        x3v = np.append(d['x3v'][0]-dph,d['x3v'])
        x3v = np.append(x3v,x3v[-1]+dph)
        var_data = np.append([d[var][-1]],d[var],axis=0)
        var_data = np.append(var_data,[var_data[0]],axis=0)
    else:
        x3v = d['x3v']
        var_data = d[var]
        
    var_interp = RegularGridInterpolator((x3v,d['x2v'],d['x1v']),var_data,bounds_error=False,method=method)
    return var_interp

def cart_to_polar(x,y,z):
    """cartesian->polar conversion (matches 0<phi<2pi convention of Athena++)
    Parameters
    x, y, z
    Returns
    r, th, phi
    """
    r = np.sqrt(x**2 + y**2 +z**2)
    th = np.arccos(z/r)
    phi = np.arctan2(y,x)
    phi = np.where(phi>=0,phi,phi+2*np.pi)
    return r,th,phi


def get_uniform_meshgrid(lim,npoints=50,center=[0,0,0]):
    """
    MM: define a uniform spacing meshgrid of cartesian points around some center
    Parameters
    ------------
    lim: float
        limit of the box, such that the cartesian box extends +/- lim relative to the center
    center: (3,) array/list of floats
        x,y,z coordinates of the center of the grid
        
    Returns
    -----------
    xx,yy,zz,rr,tt,pp: (npoints,npoints,npoints) arrays of x,y,z and r,theta,phi coordinates of the points
    
    """
    x = np.linspace(-lim,lim,npoints)+center[0]
    y = np.linspace(-lim,lim,npoints)+center[1]
    z = np.linspace(-lim,lim,npoints)+center[2]

    yy,xx,zz = np.meshgrid(y,x,z)
    r,th,ph = cart_to_polar(xx,yy,zz)

    return xx,yy,zz,r,th,ph

def get_marching_cubes_mesh(xx,yy,zz,data,level,step_size=1):
    """ MM: call marching cubes algorithm to get a surface 
    Parameters
    ----------
    xx,yy,zz: (N,N,N) arrays of x,y,z coordinates
    data: (N,N,N) array of floats in which to find surface
    level: float value to construct surface at
    step_size: down sampling of 3D data
    dx: x,y,z spacing of data
    
    Returns
    -------
    verts : (V, 3) array
        Spatial coordinates for V unique mesh vertices. Coordinate order
        matches input `volume` (M, N, P). If ``allow_degenerate`` is set to
        True, then the presence of degenerate triangles in the mesh can make
        this array have duplicate vertices.
    faces : (F, 3) array
        Define triangular faces via referencing vertex indices from ``verts``.
        This algorithm specifically outputs triangles, so each face has
        exactly three indices.
    centroids : (F,3) array
        Centroids of triangular faces
    areas : (F,) array
        areas of triangular faces
    normals : (F,3) array
        face-centered normals 
        
        
    NOTE: 
    assumes uniform spacing in x,y,z
    
    """
    center = np.mean(xx[:,0,0]), np.mean(yy[0,:,0]), np.mean(zz[0,0,:])
    lim = xx[-1,0,0]
    dx = xx[1,0,0]-xx[0,0,0]
    
    verts, faces,_,_ =  measure.marching_cubes(data,
                                               level=level,
                                               step_size=step_size,
                                               spacing=[dx,dx,dx])
    verts = verts-lim+center  # recenter the verticies
    
    centroids = mesh_get_centroids(verts,faces)
    areas     = mesh_get_areas(verts,faces)
    normals   = mesh_get_centroid_normals(verts,faces)
    
    mesh = {}
    mesh['verts'] = verts
    mesh['faces'] = faces 
    mesh['centroids'] = centroids
    mesh['areas'] = areas
    mesh['normals'] = normals
    return mesh
    #return verts,faces, centroids, areas, normals

def mesh_get_areas(verts, faces):
    """
    MM: based on skiimage routine "measure.mesh_surface_area()"
    Compute surface area, given vertices & triangular faces
    Parameters
    ----------
    verts : (V, 3) array of floats
        Array containing (x, y, z) coordinates for V unique mesh vertices.
    faces : (F, 3) array of ints
        List of length-3 lists of integers, referencing vertex coordinates as
        provided in `verts`
    Returns
    -------
    areas : (F,3) array of floats
        Surface areas of mesh triangles. Units now [coordinate units] ** 2.
    Notes
    -----
    The arguments expected by this function are the first two outputs from
    `skimage.measure.marching_cubes`. For unit correct output, ensure correct
    `spacing` was passed to `skimage.measure.marching_cubes`.
    This algorithm works properly only if the ``faces`` provided are all
    triangles.
    See Also
    --------
    skimage.measure.marching_cubes
    skimage.measure.marching_cubes_classic
    """
    # Fancy indexing to define two vector arrays from triangle vertices
    actual_verts = verts[faces]
    a = actual_verts[:, 0, :] - actual_verts[:, 1, :]
    b = actual_verts[:, 0, :] - actual_verts[:, 2, :]
    del actual_verts

    # Area of triangle in 3D = 1/2 * Euclidean norm of cross product
    #return ((np.cross(a, b) ** 2).sum(axis=1) ** 0.5).sum() / 2.
    return ((np.cross(a, b) ** 2).sum(axis=1) ** 0.5) / 2.

def mesh_get_centroids(verts,faces):
    """
    MM: from verts, faces, return the coordinates of the centroids
    
    Parameters
    ----------
    verts : (V, 3) array of floats
        Array containing (x, y, z) coordinates for V unique mesh vertices.
    faces : (F, 3) array of ints
        List of length-3 lists of integers, referencing vertex coordinates as
        provided in `verts`
    Returns
    -------
    centroids: (F, 3) array of floats
        array containing (x,y,z) coordinates of triangle centroids    
    """
    return verts[faces].sum(axis=1)/3.


def mesh_get_centroid_normals(verts,faces):
    """
    MM: from verts, faces, return the normals at the centroids
    
    Parameters
    ----------
    verts : (V, 3) array of floats
        Array containing (x, y, z) coordinates for V unique mesh vertices.
    faces : (F, 3) array of ints
        List of length-3 lists of integers, referencing vertex coordinates as
        provided in `verts`
    Returns
    -------
    normals: (F, 3) array of floats
        array containing (x,y,z) components of normal to the plane defined by the vertices of the triangle.     
    """
       # Fancy indexing to define two vector arrays from triangle vertices
    actual_verts = verts[faces]
    a = actual_verts[:, 0, :] - actual_verts[:, 1, :]
    b = actual_verts[:, 0, :] - actual_verts[:, 2, :]
    del actual_verts
    
    # cross product is perpendicular to two vectors connecting vertecies
    cp = np.cross(b, a)
    norms = (cp.T/np.linalg.norm(cp,axis=1)).T
    
    return (cp.T/np.linalg.norm(cp,axis=1)).T

def mesh_interpolate_at_xyzpoints(d,var,points):
    """
    MM: convience function to interpolate a variable to mesh points
    Parameters
    -----------
    d: athena++ data dict
    var: str variable name in, e.g. "rho"
    points: array of cartesian positions (eg vertices or centroids) (N,3) floats N x (x,y,z)
    """
    var_interp = get_interp_function(d,var)
    rp,thp,php = cart_to_polar(points[:,0],points[:,1],points[:,2])
    return var_interp( (php,thp,rp) )




def get_req_func(d,logtau):
    """
    MM: Given a dataset with variable 'tau' defined, compute the radius of the surface of logtau as a function of theta, phi
    ---------------
    d: athena++ data dict (must include tau)
    logtau: the value of logtau to interpolate radius of
    
    returns: 2d function f(theta) for radius of the star. to compute the oblate radius of the donor. 
    """
    Req = np.zeros_like(d['rho'][:,:,0])

    for k in range(len(d['x3v'])):
        for j in range(len(d['x2v'])):
            Req[k,j] = np.interp(logtau, np.flip(np.log10(d['tau'][k,j,:])), np.flip(d['gx1v'][k,j,:]) ) 
            
    ReqAzAvg = np.mean(Req,axis=0)

    return interp1d(d['x2v'],ReqAzAvg,bounds_error=False,fill_value="extrapolate")



def get_req_func_mesh(mymesh):
    r = np.linalg.norm(mymesh['verts'],axis=1)
    R = np.linalg.norm(mymesh['verts'][:,0:2],axis=1)
    z = mymesh['verts'][:,2]
    th = np.arctan(z/R) + np.pi/2

    plt.scatter(th,r)

    thf = np.linspace(0,np.pi,61)
    reqlist = []
    for i in range(len(thf)-1):
        thl = thf[i]
        thr = thf[i+1]
        sel = (th>=thl) & (th<thr)
        myr = np.mean(r[sel])
        reqlist.append([(thl+thr)/2, myr])

    reqlist = np.array(reqlist)

    plt.plot(reqlist[:,0],reqlist[:,1],'C1-')

    return interp1d(reqlist[:,0],reqlist[:,1],bounds_error=False,fill_value="extrapolate")


def h_ld(mu):
    """Eddington limb-darkening law (e.g. Pfahl2008 eq 10)"""
    return 1+1.5*mu


"""
def flux(Teff):
    return Teff

def flux_g(mesh,beta=2):
    return np.linalg.norm(mesh['centroids'],axis=1)**-beta

def Jintegral(mesh,n0):
    #Integrate J over the object surface
    ndn_all = np.dot(mesh["normals"],n0)
    ndn = np.where(ndn_all>0,ndn_all,0)
    return np.nansum(mesh["areas"]*ndn*h_ld(ndn)*flux(mesh["temp_faces"]))

def Jintegral_g(mesh,n0,beta=1):
    #Integrate J over the object surface
    ndn_all = np.dot(mesh["normals"],n0)
    ndn = np.where(ndn_all>0,ndn_all,0)
    return np.nansum(mesh["areas"]*ndn*h_ld(ndn)*flux_g(mesh,beta))
"""


def flux_xi(mesh,reqfunc,lscale=-1,M1=34.5*c.msun,Omega_env=0,alpha=1,grav_darken=True): 
    r_centroids,th_centroids,ph_centroids = cart_to_polar( mesh['centroids'][:,0], mesh['centroids'][:,1],mesh['centroids'][:,2])
    xir_o_r = mesh['xir_faces']/r_centroids
    if grav_darken:
        geff = np.abs( - 6.674e-8*M1/reqfunc(th_centroids)**2 + Omega_env**2 * reqfunc(th_centroids)*np.sin(th_centroids) )
    else:
        geff = 1.0
    flux = (1+lscale*xir_o_r) * geff
    flux = np.where(flux>0,flux,0)**alpha
    return flux

def Jintegral_xi(mesh,n0,reqfunc,lscale=-1,M1=34.5*c.msun,Omega_env=0,alpha=1,grav_darken=True):
    """Integrate J over the object surface"""
    ndn_all = np.dot(mesh["normals"],n0)
    ndn = np.where(ndn_all>0,ndn_all,0)
    return np.nansum(mesh["areas"]*ndn*h_ld(ndn)*flux_xi(mesh,reqfunc,lscale=lscale,M1=M1,Omega_env=Omega_env,alpha=alpha,grav_darken=grav_darken))



if __name__=='__main__':
    
    ###### SIMULATION PARAMS   #########
    parser = argparse.ArgumentParser(description='Read input/output directories')
    parser.add_argument("--base_dir", help="data directory (should end with / )")
    parser.add_argument("--level",help="refinement level at which to read data")
    parser.add_argument("--nxyz",help="number of mesh zones in xyz remesh",default=200,type=int)
    parser.add_argument("--logrho",help='surface of constant tau to extract',default=0.0,type=float,nargs='+')
    parser.add_argument("--filestart",help='int file number to start at',default=0,type=int)
    parser.add_argument("--filestop",help='int file number to stop at',default=-1,type=int)
    parser.add_argument("--fileskip",help='int file number to skip',default=1,type=int)
    parser.add_argument("--lim",help="xyz limits of surface search",default=3e12,type=float)

    args = parser.parse_args()
    base_dir=args.base_dir
    if args.level:
        mylevel = int(args.level)
    nxyz=args.nxyz
    mylim=args.lim
    logrho = args.logrho
    filestart=args.filestart
    filestop =args.filestop
    fileskip =args.fileskip

    file_list = sorted(glob(base_dir+"HSE.out1.[0-9][0-9][0-9][0-9][0-9].athdf"))
    file_list = file_list[filestart:filestop:fileskip]
    print( file_list )
    print("slicing at logrho:", logrho)
    ####################################

    orb = ou.read_trackfile(base_dir+"pm_trackfile.dat")

    xx,yy,zz,r,th,ph = get_uniform_meshgrid(mylim,npoints=nxyz,center=[0,0,0])
        
    def read_file_save_mesh(fn):
        d=ou.read_data(fn,orb,x1_max=2*mylim,level=mylevel,gamma=1.35)
        # define interpolating functions & get data
        d['temp'] = 0.61*c.mp/c.kB * d['press']/d['rho']
        #M1 = np.sum(d['rho']*d['r0']*d['dvol']) + orb['m1'][0]
        #d['g'] = (6.674e-8*M1/d['gx1v']**2)
        #Hp = d['press']/(d['g']*d['rho'])
        #d["tau"] = d['rho']*np.where(d['r0']>0.9,d['r0'],0) *kappa(d['rho'],d['temp'])*Hp

        d['rhor0'] = d['rho']*d['r0']
        rho_interp = get_interp_function(d,'rhor0')
        rho_vals = rho_interp((ph,th,r))


        for lr in logrho:
            # get a surface mesh object (constant tau)
            m = get_marching_cubes_mesh(xx,yy,zz,np.log10(rho_vals),lr,step_size=1)
            
            # interpolate other values to the vertices of the mesh
            m["r0_faces"]    = mesh_interpolate_at_xyzpoints(d,"r0",m["centroids"])
            m["temp_faces"]  = mesh_interpolate_at_xyzpoints(d,"temp",m["centroids"])
            m["rho_faces"]   = mesh_interpolate_at_xyzpoints(d,"rho",m["centroids"])
            m["press_faces"] = mesh_interpolate_at_xyzpoints(d,"press",m["centroids"])
            m["r0_verts"]    = mesh_interpolate_at_xyzpoints(d,"r0",m["verts"])
            m["temp_verts"]  = mesh_interpolate_at_xyzpoints(d,"temp",m["verts"])
            m["rho_verts"]   = mesh_interpolate_at_xyzpoints(d,"rho",m["verts"])
            m["press_verts"] = mesh_interpolate_at_xyzpoints(d,"press",m["verts"])
            m['vel1_verts']  = mesh_interpolate_at_xyzpoints(d,"vel1",m["verts"])
            m['vel2_verts']  = mesh_interpolate_at_xyzpoints(d,"vel2",m["verts"])
            m['vel3_verts']  = mesh_interpolate_at_xyzpoints(d,"vel3",m["verts"])
            
            # metadata
            m['time'] = d['Time']
        
            # save the data
            dd.io.save(base_dir+"mesh_"+fn[-11:-6]+"_rho1e"+str(lr)+".h5",m)
        return fn

    #ml = []
    #for fn in file_list:
    #    ml.append( read_file_return_mesh(fn) )
    
    ntasks = os.environ['SLURM_NTASKS']
    with mp.Pool(int(ntasks)) as pool:
        ml = pool.map(read_file_save_mesh, file_list)


    print("completed processing: ", ml)

