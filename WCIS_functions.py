import numpy as np
import pandas as pd
import mdtraj as md
from scipy.spatial import cKDTree
from skimage import measure
import sys
from time import time
import itertools

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
    sel = ' or '.join(['name '+m for mol in molecules for m in mol.split(' ') if m[0]!='H'])
    traj = md.load_xtc(filename,top=top,
            atom_indices=top.top.select(sel))[nskip:]
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
                                pair=atoms, theta=np.zeros(theta.size), all=np.empty(0), ndang=np.empty(0)), 
                              'lower': dict(cosine=np.zeros(z.size), conc=np.zeros(z.size), 
                                pair=atoms, theta=np.zeros(theta.size), all=np.empty(0), ndang=np.empty(0))}
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
             "z": z,
             "sel": sel,
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
        pos = frame.atom_slice(frame.top.select('name '+atom)).xyz[0]
        _, ind = tree.query(pos, k=1)
        dist = (tree.data[ind,2] - pos[:,2])*sign
        if type(value[label]) == dict:
            hist, _ = np.histogram(dist,bins=edges,density=False)
            value[label]['conc'] += hist * toM
            if value[label]['theta'].ndim > 1:
                cosine(frame,dist,atom,normals[ind,:],value[label],sign,edges,toM,layers)
            else:
                selOSN = frame.top.select('name O or name S1 or name N3')
                posOSN = frame.atom_slice(selOSN).xyz[0]
                _, indOSN = tree.query(posOSN, k=1)
                distOSN = (tree.data[indOSN,2] - posOSN[:,2])*sign
                mask = np.logical_and(dist>-.3,dist<.3)
                idx = np.asarray(frame.top.select('name O'))[mask]
                ind_surf = ind[mask]
                mask = np.logical_and(distOSN>-.3,distOSN<.7)
                idxOSN = np.asarray(selOSN)[mask]
                dangling(frame,value[label],sign,dist,distOSN,idx,idxOSN,normals[ind_surf,:])     
        else:
            hist, _ = np.histogram(dist,bins=edges,density=False)
            value[label] += hist * toM

def cosine(frame,dist,atom,normals,dictionary,sign,edges,toM,layers):
    selection_string = 'name '+dictionary['pair'][0]+' or name '+dictionary['pair'][1]
    pair = np.array(frame.top.select(selection_string)).reshape(-1,2)
    vec = md.compute_displacements(frame,pair).reshape(-1,3)
    cosine = np.einsum('ij,ij->i',vec,normals) / np.linalg.norm(vec,axis=1)
    if atom=='C2':
        cosine = -cosine
    hist, _ = np.histogram(dist,bins=edges,weights=cosine,density=False)
    dictionary['cosine'] += hist * toM
    angle = np.arccos(np.clip(cosine,-1,1))/np.pi*180
    for i in range(len(layers)-1):
        mask = np.logical_and(dist>layers[i],dist<layers[i+1])
        hist, _ = np.histogram(angle[mask],bins=np.arange(0,181,2),density=False)
        dictionary['theta'][i] += hist

def dangling(frame,dictionary,sign,dist,distOSN,idx,idxOSN,normals):
    pair_ox = np.asarray(list(itertools.product(idx, idxOSN)))
    # remove pairs of same index
    pair_ox = pair_ox[np.std(pair_ox,axis=1)!=0]
    dist_ox = md.compute_distances(frame,pair_ox).reshape(idx.size,-1)
    h1o = np.c_[idx,idx+1]
    h2o = np.c_[idx,idx+2]
    # compute H-O vectors
    vec_h1o = md.compute_displacements(frame,h1o).reshape(-1,3)
    vec_h2o = md.compute_displacements(frame,h2o).reshape(-1,3)
    cosine_h1o = np.einsum('ij,ij->i',vec_h1o,normals) / np.linalg.norm(vec_h1o,axis=1)
    cosine_h2o = np.einsum('ij,ij->i',vec_h2o,normals) / np.linalg.norm(vec_h2o,axis=1)
    angle_h1o = np.arccos(np.clip(cosine_h1o,-1,1))/np.pi*180
    angle_h2o = np.arccos(np.clip(cosine_h2o,-1,1))/np.pi*180
    h1ox = np.c_[np.asarray(pair_ox)[:,0]+1,np.asarray(pair_ox)[:,0],np.asarray(pair_ox)[:,1]]
    h2ox = np.c_[np.asarray(pair_ox)[:,0]+2,np.asarray(pair_ox)[:,0],np.asarray(pair_ox)[:,1]]
    # compute H-O...O angles
    angle_h1ox = md.compute_angles(frame,h1ox).reshape(idx.size,-1)/np.pi*180
    angle_h2ox = md.compute_angles(frame,h2ox).reshape(idx.size,-1)/np.pi*180
    # selection of dangling OH based on R-beta definition (DOI: 10.1021/acs.jctc.7b00566)
    Rc = dist_ox<=.35
    beta1 = angle_h1ox<=50
    beta2 = angle_h2ox<=50
    angles_ho = np.append(angle_h1o[np.sum(beta1*Rc,axis=1)==0],angle_h2o[np.sum(beta2*Rc,axis=1)==0])
    dictionary['ndang'] = np.append(dictionary['ndang'],angles_ho.size)
    dictionary['all'] = np.append(dictionary['all'],idx.size)
    hist, _ = np.histogram(angles_ho,bins=np.arange(0,181,2),density=False)
    dictionary['theta'] += hist

def find_isosurfaces(frame,params):
    radius,factor,scale,sel = params['radius'],params['factor'],params['scale'],params['sel']
    level = params['level']
    params['box'] = frame.unitcell_lengths[0]
    params['grid'],params['spacing'],params['grid_shape'] = make_grid(params['box'],params['mesh'])
    grid,spacing,grid_shape,box = params['grid'],params['spacing'],params['grid_shape'],params['box']
    
    pos = frame.atom_slice(atom_indices=frame.top.select(sel)).xyz[0]
            
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