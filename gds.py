import numpy as np
import pandas as pd
import mdtraj as md
from scipy.optimize import curve_fit
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

def initialize(filename,top,bw,molecules,nskip,layers):
    traj = md.load_xtc(filename,top=top,atom_indices=top.top.select('name O'))[nskip:]
    box = traj.unitcell_lengths.mean(0)
    com = md.compute_center_of_mass(traj)
    z_com = np.repeat(com[:,2],traj.n_atoms)
    z_com = z_com.reshape(traj.n_frames,-1)
    z_pos = traj.xyz[:,:,2] - z_com
    hist, edges = np.histogram(z_pos,bins=np.arange(-box[2]*.5,box[2]*.5,bw))
    z = edges[:-1]+(edges[1]-edges[0])/2.
    hist = hist/(bw*box[0]*box[1]*6.022*0.1*traj.n_frames)
    
    z_min = np.floor(z[hist>hist.max()/2].min()).astype(int)
    z_max = np.ceil(z[hist>hist.max()/2].max()).astype(int)
    z_mid = (z_max+z_min)*.5
    rho0 = hist[np.logical_and(z>z_min+2,z<z_max-2)].mean()
        
    hyptan = lambda x,gds,d : rho0/2*(1-np.tanh((x-gds)/d))

    popt,pcov = curve_fit(hyptan,z[z>z_mid],hist[z>z_mid],p0=[z_max,.1])
    z_gds1 = popt[0]
    thickness1 = popt[1]

    popt,pcov = curve_fit(hyptan,-z[z<z_mid],hist[z<z_mid],p0=[-z_min,.1])
    z_gds0 = -popt[0]
    thickness2 = popt[1]
        
    level = 0.5*hist[np.logical_and(z>z_min+2,z<z_max-2)].mean()
    edges = np.arange(-1,z_max-z_min+2,bw*4)
    z = edges[:-1]+(edges[1]-edges[0])/2.
    thetaedges = np.arange(0,181,2)
    theta = thetaedges[:-1]+(thetaedges[1]-thetaedges[0])/2.
    data = {}
    for molecule in molecules:
        atoms = molecule.split(' ')
        print(atoms)
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
    params = {"box": box,
             "bw": bw*4,
             "n_frames": traj.n_frames+nskip,
             "edges": edges,
             "theta": theta,
             "z": z,
             "level": level,
             "layers": layers,
             "z_gds": [z_gds0,z_gds1],
             "thickness": [thickness1,thickness2],
             "data": data}
    del traj
    return pd.Series(params)

def calc_profiles(frame,params,label,sign,layers):
    box, edges, data, bw = params['box'], params['edges'], params['data'], params['bw']
    z_gds = params['z_gds'][(sign+1)//2]
    toM = 1. / (bw*box[0]*box[1]*6.022*0.1)
    for atom, value in data.items():
        oxygens = frame.top.select('name O')
        z_com = md.compute_center_of_mass(frame.atom_slice(oxygens))[:,2]
        z_pos = frame.atom_slice(frame.top.select('name '+atom)).xyz[0,:,2] - z_com
        dist = (z_gds - z_pos)*sign
        if type(value[label]) == dict:
            hist, _ = np.histogram(dist,bins=edges,density=False)    
            value[label]['conc'] += hist * toM
            profile_cosine(frame,dist,atom,value[label],sign,edges,toM,layers)
            if value[label]['theta'].ndim == 1:
                selOSN = frame.top.select('name O or name S1 or name N3')
                z_pos = frame.atom_slice(selOSN).xyz[0,:,2] - z_com
                distOSN = (z_gds - z_pos)*sign
                mask = np.logical_and(dist>-.3,dist<.3)
                idx = np.asarray(oxygens)[mask]
                mask = np.logical_and(distOSN>-.3,distOSN<.7)
                idxOSN = np.asarray(selOSN)[mask]
                dangling(frame,value[label],sign,dist,distOSN,idx,idxOSN)
        else:
            hist, _ = np.histogram(dist,bins=edges,density=False)
            value[label] += hist * toM    

def profile_cosine(frame,dist,atom,dictionary,sign,edges,toM,layers):
    selection_string = 'name '+dictionary['pair'][0]+' or name '+dictionary['pair'][1]
    pair = np.array(frame.top.select(selection_string)).reshape(-1,2)
    vec = md.compute_displacements(frame,pair).reshape(-1,3)
    cosine = vec[:,2] / np.linalg.norm(vec,axis=1) * sign
    if atom=='C2':
        cosine = -cosine
    hist, _ = np.histogram(dist,bins=edges,weights=cosine,density=False)
    dictionary['cosine'] += hist * toM
    if dictionary['theta'].ndim > 1:
        angle = np.arccos(np.clip(cosine,-1,1))/np.pi*180
        for i in range(len(layers)-1):
            mask = np.logical_and(dist>layers[i],dist<layers[i+1])
            hist, _ = np.histogram(angle[mask],bins=np.arange(0,181,2),density=False)
            dictionary['theta'][i] += hist

def dangling(frame,dictionary,sign,dist,distOSN,idx,idxOSN):
    pair_ox = np.asarray(list(itertools.product(idx, idxOSN)))
    # remove pairs of same index
    pair_ox = pair_ox[np.std(pair_ox,axis=1)!=0]
    dist_ox = md.compute_distances(frame,pair_ox).reshape(idx.size,-1)
    h1o = np.c_[idx,idx+1]
    h2o = np.c_[idx,idx+2]
    # compute H-O vectors
    vec_h1o = md.compute_displacements(frame,h1o).reshape(-1,3)
    vec_h2o = md.compute_displacements(frame,h2o).reshape(-1,3)
    cosine_h1o = vec_h1o[:,2] / np.linalg.norm(vec_h1o,axis=1) * sign
    cosine_h2o = vec_h2o[:,2] / np.linalg.norm(vec_h2o,axis=1) * sign
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

t1 = time()
print('molecules',args.molecules)
print('layers',args.layers)
top = md.load(args.top)
params = initialize(args.traj,top,0.005,args.molecules,args.begin,args.layers)
cnt = 0
for i in range(args.begin,params['n_frames'],args.every):
    frame = md.load_frame(args.traj, i, top=top)
    calc_profiles(frame,params,label='upper',sign=1,layers=args.layers)
    calc_profiles(frame,params,label='lower',sign=-1,layers=args.layers)
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
