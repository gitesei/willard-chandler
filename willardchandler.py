import numpy as np
import pandas as pd
import mdtraj as md
from scipy.spatial import cKDTree
from skimage import measure
import sys
from time import time
import itertools
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('-s',dest='top')
parser.add_argument('-t',dest='traj')
parser.add_argument('-b',dest='begin',type=int)
parser.add_argument('-e',dest='every',type=int)
parser.add_argument('-m',dest='molecules',nargs='+')
parser.add_argument('-l',dest='layers',nargs='+',type=float)
parser.add_argument('-o',dest='out')
args = parser.parse_args()

def make_grid(box,mesh):
    ngrid = np.ceil(box / mesh).astype(int)
    grid_shape = tuple(ngrid.astype(int))
    spacing = box / ngrid
    xyz = []
    for i in range(3):
        xyz.append(np.linspace(0., box[i]-box[i]/ngrid[i], ngrid[i]))
    x, y, z = np.meshgrid(xyz[0], xyz[1], xyz[2], indexing='ij')
    grid = np.c_[x.reshape(-1, 1), y.reshape(-1, 1), z.reshape(-1, 1)]
    return grid, spacing, grid_shape

def initialize(filename,top,mesh,alpha,molecules,nskip,layers):
    traj = md.load_xtc(filename,top=top,
            atom_indices=top.top.select('all and not name H'))[nskip:]
    box = traj.unitcell_lengths.mean(0)    
    hist, edges = np.histogram(traj.xyz[:,:,2],bins=np.arange(0,box[2],1).astype(int))  
    z = edges[:-1]+(edges[1]-edges[0])/2.
    hist = hist/(box[0]*box[1]*traj.n_frames)
    z_min = np.floor(z[hist>hist.max()/2].min()).astype(int)
    z_max = np.ceil(z[hist>hist.max()/2].max()).astype(int)
    level = 0.5*hist[np.logical_and(z>z_min+2,z<z_max-2)].mean()
    edges = np.arange(-1,z_max-z_min+2,mesh*0.1)
    z = edges[:-1]+(edges[1]-edges[0])/2.
    thetaedges = np.arange(0,181,2)
    theta = thetaedges[:-1]+(thetaedges[1]-thetaedges[0])/2.
    grid, spacing, grid_shape = make_grid(box,mesh)
    data = {}
    for molecule in molecules:
        atoms = molecule.split(' ')
        if len(atoms)==1:
            data[atoms[0]] = {'upper': np.zeros(z.size), 'lower': np.zeros(z.size)}
        elif len(atoms)==3:
            for atom in atoms[::2]:
                data[atom] = {'upper': np.zeros(z.size), 'lower': np.zeros(z.size)}
            data[atoms[1]] = {'upper': dict(cosine=np.zeros(z.size), conc=np.zeros(z.size), 
                                pair=atoms[::2], theta=np.zeros((len(layers)-1,theta.size))), 
                              'lower': dict(cosine=np.zeros(z.size), conc=np.zeros(z.size), 
                                pair=atoms[::2], theta=np.zeros((len(layers)-1,theta.size)))}
        elif len(atoms)==2:
            data[atoms[0]] = {'upper': np.zeros(z.size), 'lower': np.zeros(z.size)}
            data[atoms[1]] = {'upper': dict(cosine=np.zeros(z.size), conc=np.zeros(z.size), 
                                pair=atoms, theta=np.zeros((len(layers)-1,theta.size))), 
                              'lower': dict(cosine=np.zeros(z.size), conc=np.zeros(z.size), 
                                pair=atoms, theta=np.zeros((len(layers)-1,theta.size)))}
    params = {"radius": 3*alpha,
             "scale": 2.*alpha**2,
             "factor": np.power(2.*np.pi*alpha**2,-3/2.),
             "box": box,
             "mesh": mesh,
             "bw": mesh*0.1,
             "n_frames": traj.n_frames+nskip,
             "spacing": spacing,
             "grid": grid,
             "grid_shape": grid_shape,
             "edges": edges,
             "theta": theta,
             "thetaedges": thetaedges,
             "z": z,
             "surface_area": np.empty(0),
             "surface_zstd": np.empty(0),
             "level": level,
             "layers": layers,
             "data": data}
    del traj
    return pd.Series(params)

def calc_profiles(frame,surface,params,label,sign,layers):
    box, edges, data, bw = params['box'], params['edges'], params['data'], params['bw']
    toM = 1. / (bw*box[0]*box[1]*6.022*0.1)
    verts, normals = surface
    tree = cKDTree(verts, boxsize=box)
    for atom, value in data.items():
        if type(value[label]) == dict:
            profile_cosine(frame,atom,surface,tree,value[label],sign,edges,toM,layers)
        else:
            profile(frame,atom,tree,value[label],sign,edges,toM)

def profile(frame,atom,tree,array,sign,edges,toM):
    pos = frame.atom_slice(frame.top.select('name '+atom)).xyz[0]
    _, ind = tree.query(pos, k=1)
    dist = (tree.data[ind,2] - pos[:,2])*sign
    hist, _ = np.histogram(dist,bins=edges,density=False)
    array += hist * toM

def profile_cosine(frame,atom,surface,tree,dictionary,sign,edges,toM,layers):
    verts, normals = surface
    pos = frame.atom_slice(frame.top.select('name '+atom)).xyz[0]
    _, ind = tree.query(pos, k=1)
    dist = (tree.data[ind,2] - pos[:,2])*sign
    hist, _ = np.histogram(dist,bins=edges,density=False)
    dictionary['conc'] += hist * toM
    selection_string = 'name '+dictionary['pair'][0]+' or name '+dictionary['pair'][1]
    pair = np.array(frame.top.select(selection_string)).reshape(-1,2)
    vec = md.compute_displacements(frame,pair).reshape(-1,3)
    cosine = np.einsum('ij,ij->i',vec,normals[ind,:]) / np.linalg.norm(vec,axis=1)
    # if atom=='C2':
    #     cosine = -cosine
    hist, _ = np.histogram(dist,bins=edges,weights=cosine,density=False)
    dictionary['cosine'] += hist * toM
    angle = np.arccos(np.clip(cosine,-1,1))/np.pi*180
    thetaedges = np.arange(0,181,2)
    theta = thetaedges[:-1]+(thetaedges[1]-thetaedges[0])/2.
    for i in range(len(layers)-1):
        mask = np.logical_and(dist>layers[i],dist<layers[i+1])
        hist, _ = np.histogram(angle[mask],bins=thetaedges,density=False)
        dictionary['theta'][i] += hist

def find_isosurfaces(frame,params):
    radius,factor,scale = params['radius'],params['factor'],params['scale']
    level = params['level']
    params['box'] = frame.unitcell_lengths[0]
    params['grid'],params['spacing'],params['grid_shape'] = make_grid(params['box'],params['mesh'])
    grid,spacing,grid_shape,box = params['grid'],params['spacing'],params['grid_shape'],params['box']
    
    pos = frame.atom_slice(atom_indices=frame.top.select('name O or resname SCN or name NA')).xyz[0]
            
    tree = cKDTree(grid,boxsize=box)
                
    # indeces of grid points within a radial distance from the particles
    indlist = tree.query_ball_point(pos, radius)

    # unwrapped list of lists
    indarray = np.asarray(list(itertools.chain.from_iterable(indlist)),dtype=int)

    # lenghts of the sublists in indlist
    lenarray = np.asarray([len(ind) for ind in indlist],dtype=int)

    # vector distance between particles and grid points
    dr = grid[indarray,:] - np.repeat(pos, lenarray, axis=0)

    # periodic boundary conditions in xy-plane
    cond = np.where(np.abs(dr) > box / 2.)
    dr[cond] -= np.sign(dr[cond])*box[cond[1]]
    # coarse grained density field
    dens = factor*np.exp( - np.linalg.norm(dr,ord=2,axis=1)**2 / scale )
    # densities at the same grid point are summed up
    field = pd.DataFrame(data={'index': indarray, 'dens': dens}).groupby('index').sum()
    # grid points with zero density are included in the dataframe
    new_index = pd.Index(range(grid.shape[0]), name="index")
    field = field.reindex(new_index,fill_value=0).values.reshape(grid_shape)
    
    verts, faces, normals, values = measure.marching_cubes_lewiner(field, level, spacing=tuple(spacing))
    verts[:,:2] = verts[:,:2] + spacing[:2] / 2.
    cond_upper = verts[:,2]>box[2]/2
    upper = verts[cond_upper], normals[cond_upper]
    lower = verts[~cond_upper], normals[~cond_upper]
    params['surface_area'] = np.append(params['surface_area'],measure.mesh_surface_area(verts,faces)*0.5)
    params['surface_zstd'] = np.append(params['surface_zstd'],upper[0][:,2].std())
    params['surface_zstd'] = np.append(params['surface_zstd'],lower[0][:,2].std())
    return upper, lower

t1 = time()
print('molecules',args.molecules)
print('layers',args.layers)
top = md.load(args.top)
params = initialize(args.traj,top,0.2,0.24,args.molecules,args.begin,args.layers)
cnt = 0
for i in range(args.begin,params['n_frames'],args.every):
    frame = md.load_frame(args.traj, i, top=top)
    upper, lower = find_isosurfaces(frame,params)
    calc_profiles(frame,upper,params,label='upper',sign=1,layers=args.layers)
    calc_profiles(frame,lower,params,label='lower',sign=-1,layers=args.layers)
    cnt += 1
for molecule in args.molecules:
    atoms = molecule.split(' ')
    for atom in atoms:
        if type(params['data'][atom]['upper']) == dict:
            params['data'][atom]['upper']['conc'] /= cnt 
            params['data'][atom]['lower']['conc'] /= cnt 
            params['data'][atom]['upper']['cosine'] /= cnt 
            params['data'][atom]['lower']['cosine'] /= cnt 
        else:
            params['data'][atom]['upper'] /= cnt 
            params['data'][atom]['lower'] /= cnt
pd.to_pickle(params,args.out)
t2 = time()
print((t2-t1)/3600)
