'''
Created on March 20, 2020
@author: Ruhar 
'''
#! /usr/local/bin/env python3

#IMPORTING ALL THE MSMBUILDER REQUIRED MODULES 
import numpy as np                                                                                                                                                                                                                                                                                                                                      
import sys                                                                                                                                                                                          
import math                                                                                                                                              
import csv                                                                                                                                                
import multiprocessing as multiproc                                                                                                                        
import time                                                                                        
from argparse import ArgumentParser                                                                                                      
from msmbuilder.dataset import dataset                                 
from msmbuilder.featurizer import AlphaAngleFeaturizer                 
from msmbuilder.decomposition import tICA    
from msmbuilder.cluster import MiniBatchKMeans      
from msmbuilder.msm import ContinuousTimeMSM        
from msmbuilder.utils import verbosedump,verboseload
from msmbuilder.utils.nearest import KDTree                                                                                             
import os,glob                                                                                                                          
import numpy as np                                                                                                  
import mdtraj as md                                                                                                 
import pickle                                                                                                           
import seaborn as sns                                                                                                                                                 
import numpy as np                                                                                                                                                      
from matplotlib import pyplot as plt                                                                                    
from msmbuilder.featurizer import DihedralFeaturizer
from scipy.constants import Avogadro, Boltzmann, calorie_th
import matplotlib
import matplotlib.patches as mpatches
from seaborn.distributions import (_scipy_univariate_kde, _scipy_bivariate_kde)

###READ TRAJECTORY FILES AND PREPROCESSING the data
ds=dataset("*.nc", topology="s.pdb")                                                                                     
feat = DihedralFeaturizer(types=['phi', 'psi'])                                                                                                                       
ds_alpha=ds.fit_transform_with(feat, "dihed/",fmt='dir-npy')                                                                                                           
ds_alpha = dataset("./dihed/")                                                                                         
print(len(ds_alpha),len(ds))                                                                                           
print(ds[0].xyz.shape)                                                              
ds_alpha = dataset("./dihed/")                                                                                          
tica_mdl = tICA(lag_time=10,n_components=2)                                                                                                                          
tica_features = ds_alpha.fit_transform_with(tica_mdl, out_ds = 'tica')                                                                                                 
tica_features = dataset("./tica/")                                                                                     
kmeans_mdl = MiniBatchKMeans(10)                                                                                       
assignments = tica_features.fit_transform_with(kmeans_mdl, out_ds='assignments/')                                       
assignments = dataset("assignments/")  
import msmexplorer as msme                                                                                                                                                                          
from matplotlib import pyplot as plt        
tica_trajs=dataset("./tica/")

###PREPROCESSING THE TICA DATA
j=0
tica=[]
for k in range(6):
    f1=[]
    for i in range(len(tica_trajs)):
        f=list(tica_trajs[i][j:j+800])
        f1.append(f)
    f1=np.array(f1)
    da=np.concatenate(f1, axis=0)
    loop = range(len(da))
    str_list=list()
    da_int=da.astype(int)
    for x in loop:
        a=str(da_int[x,:])
        str_list.append(a)
    unique=list(set(str_list))
    index=list()
    number_unique=len(unique)
    print("unique strings: ", number_unique, "of total strings:", len(str_list))
    tic=[]
    for i in range(len(da)):
        tic.append(da[i])
    tica.append(tic)
    j=j+800
###CONCATENATED THE ALL TICA VALUES
tica=np.array(tica)
tica=np.concatenate(tica, axis=0)

###PLOTTING COMMON FREE ENERGY LANDSCAPE FOR ALL PH AND EXTRACTING REPRESENTATIVE CONFORMATIONS FOR MACROSTATES
temperature=300
assignments = dataset("assignments/")
from msmbuilder.msm import ContinuousTimeMSM  
msm_mdl = ContinuousTimeMSM(lag_time=10,ergodic_cutoff=1/10) 
assignments.fit_with(msm_mdl)  
msm_mdl.percent_retained_

### CALCULATING 10 MOST POPULATED CONFORMATIONS 
f=open("msm_mdl.pkl",'wb')
pickle.dump({k: msm_mdl.__dict__[k] for k in ['lag_time','n_timescales','ergodic_cutoff','verbose','sliding_window','guess','theta_','ratemat_','transmat_','countsmat_','n_states_',
'mapping_','populations_','information_','loglikelihoods_','eigenvalues_','left_eigenvectors_','right_eigenvectors_','percent_retained_'] }, f)
f.close()
q=msm_mdl.__dict__.get('populations_')
ind=np.argpartition(q, -10)[-10:]
p=ind[np.argsort(q[ind])]
print("The most populated conformations:", q[p])

### EXTRACTING REPRESENTATIVE CONFORMATIONS OF MACROSTATE
ktree=KDTree(tica_features)
trj_ds = dataset("*.nc", topology="s.pdb")
trj_list = []
for cind in p:
    pt = kmeans_mdl.cluster_centers_[cind]
    _,(t,f) = ktree.query(pt)
    print("The trajectory and frame number of representative macrostate:",t,f)
    trj_list.append(trj_ds[t][f])

#CALCULATING COMMON FREE ENERGY LANDSCAPE FOR ALL PH
pi_0=msm_mdl.__dict__.get('populations_')[np.concatenate(assignments, axis=0)] 
data = np.concatenate(tica_trajs, axis=0) 
clip = [(-np.inf, np.inf), (-np.inf, np.inf)]
levels = np.linspace(0,2,21)
def _thermo_transform(Z, temperature):  
  return - THERMO_CONSTANT * temperature * np.log(Z)
THERMO_CONSTANT = 10**-3 * Boltzmann * Avogadro / calorie_th
bw='scott'
gridsize=100
cut=3
obs=(0,1) 
prune = tica[:, obs]  
X, Y, Z = _scipy_bivariate_kde(tica[:, 0], tica[:, 1], bw, gridsize, cut, clip)
Z = _thermo_transform(Z, temperature)
n_samples=100000
bw='scott'
gridsize=100
cut=3 
vmax = np.percentile(Z, 50)
vmin = -1E-12
xticks=np.arange(0.95, -4.5, -0.58)
yticks=np.arange(4.2, -4.9, -0.93)
ax=plt.axes()
ax.contourf(X, Y, Z - Z.min(), cmap=plt.get_cmap("hsv"),levels=np.linspace(vmin, vmax, 10),alpha=0.5, vmin=vmin, vmax=vmax) 
ax.set(xlim=(-4.5, 0.95), ylim=(-4.9, 4.2), xlabel='tIC1', ylabel='tIC2')  
ax.set_xticks(xticks) 
ax.set_yticks(yticks)
ax.grid(True)
ax.grid(color='k', linestyle='-') 
cf=ax.contourf(X, Y, Z - Z.min(), cmap=plt.get_cmap("hsv"), levels=np.linspace(vmin, vmax, 10), alpha=0.5, vmin=vmin, vmax=vmax)
cbar_kwargs={'format': '%.1f', 'label': 'Free energy (kcal/mol)'}
clabel_kwargs={'fmt': '%.1f'}
clabel=True 
cbar=True 
shade=True
if cbar:
    if not cbar_kwargs:
        cbar_kwargs = {}
    if shade:
        mappable=cf
    else:
        mappable=cs
    plt.colorbar(mappable, **cbar_kwargs)
    
if clabel:                                                                                                                                                                                                
     if not cbar_kwargs:                                                                                                                                 
         clabel_kwargs={}
pos = dict(zip(range(kmeans_mdl.n_clusters), kmeans_mdl.cluster_centers_))
msme.plot_msm_network(msm_mdl, pos=pos, node_color='pomegranate', edge_color='carbon')
plt.show()

###QUANTIFICATION OF CONFORMATIONS
from msmbuilder.dataset import dataset
def row_based_idx(num_rows, num_cols, idx):
     return np.arange(1, num_rows*num_cols + 1).reshape((num_rows, num_cols))[::-1].transpose().flatten()[idx-1]

assignments=dataset("assignments/")
tica_features = dataset("./tica/")
hist3, newedgesX3, newedgesY3 = np.histogram2d(tica[:,0][0:24000], tica[:,1][0:24000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])                                                                                                    
hist4, newedgesX4, newedgesY4 = np.histogram2d(tica[:,0][24001:48000], tica[:,1][24001:48000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])                                               
hist5, newedgesX5, newedgesY5 = np.histogram2d(tica[:,0][48001:72000], tica[:,1][48001:72000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])                                                                                            
hist6, newedgesX6, newedgesY6 = np.histogram2d(tica[:,0][72001:96000], tica[:,1][72001:96000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])                                                 
hist7, newedgesX7, newedgesY7 = np.histogram2d(tica[:,0][96001:120000], tica[:,1][96001:120000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])                                              
hist8, newedgesX8, newedgesY8 = np.histogram2d(tica[:,0][120001:144000], tica[:,1][120001:144000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])                                          
      
n1=np.concatenate(hist3, axis=0)                                                                   
n2=np.concatenate(hist4, axis=0)                                                                                                         
n3=np.concatenate(hist5, axis=0)                                       
n4=np.concatenate(hist6, axis=0)                                       
n5=np.concatenate(hist7, axis=0)             
n6=np.concatenate(hist8, axis=0)   

#PLOTTING BINS VALUES FROM COMMON FREE ENERGY LANDSCAPE
fig=plt.figure()                           
fig.subplots_adjust(hspace=0.7, wspace=0.7) 
fig.text(0.5, 0.04, 'Temperature', ha='center')
fig.text(0.08, 0.5, 'Population', va='center', rotation='vertical') 
tr=[]                                      
for i in range(100):                                                                                     
    tr=[n1[i], n2[i], n3[i], n4[i], n5[i], n6[i]]
    row_based_plot_idx = row_based_idx(10, 10, i+1)
    ax=fig.add_subplot(10, 10, row_based_plot_idx)
    ax.plot(tr, color="blue")  
    ax.xaxis.set_ticks(np.arange(0,6,1))
    matplotlib.rc('xtick', labelsize=8)
    matplotlib.rc('ytick', labelsize=8)
plt.subplots_adjust(right=0.95)  
plt.show()

###pH BASED QUANTIFICATION OF CONFORMATION IN BASIN
j=0
tica=[]                                                                                                                                                
for k in range(6):       
    f1=[]                                                                                                
    for i in range(int(len(tica_trajs)/3)):                                                                                              
        f=list(tica_trajs[3*i][j:j+800])                      
        f1.append(f)                                                   
        print(3*i)                           
        j=j                                                                                                                                             
    f1=np.array(f1)                                                                                                                     
    da = np.concatenate(f1, axis=0)                                                                                 
    loop = range(len(da))                                                                                           
    str_list=list()                                                                                                                                                                                  
  
    da_int=da.astype(int)                                                                                               
    for x in loop:                                                                                                                                                                                   
                                                                                                                     
        a=str(da_int[x,:])                                                                                                                                            
        str_list.append(a)                                                                                                                                             
    unique=list(set(str_list))                                                                                         
    index=list()                                                                                                       
    number_unique = len(unique)                                                     
    print ("unique strings:  ", number_unique, "of total strings:", len(str_list))                                                                                                                   
                               
    tic=[]                                                                                                             
    for i in range(len(da)):                                                                                            
        tic.append(da[i])                       
    tica.append(tic)
    j=j+800

###CONCATENATION
tica=np.array(tica)
tica=np.concatenate(tica, axis=0)
print(np.shape(tica))

###QUANTIFICATION OF CONFORMATIONS
histt2, newedgeX, newedgeY = np.histogram2d(tica[:,0], tica[:,1], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histt3, newedgeX3, newedgeY3 = np.histogram2d(tica[:,0][0:8000], tica[:,1][0:8000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histt4, newedgeX4, newedgeY4 = np.histogram2d(tica[:,0][8001:16000], tica[:,1][8001:16000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histt5, newedgeX5, newedgeY5 = np.histogram2d(tica[:,0][16001:24000], tica[:,1][16001:24000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histt6, newedgeX6, newedgeY6 = np.histogram2d(tica[:,0][24001:32000], tica[:,1][24001:32000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histt7, newedgeX7, newedgeY7 = np.histogram2d(tica[:,0][32001:40000], tica[:,1][32001:40000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histt8, newedgeX8, newedgeY8 = np.histogram2d(tica[:,0][40001:48000], tica[:,1][40001:48000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])

n1=np.concatenate(histt3, axis=0)
n2=np.concatenate(histt4, axis=0) 
n3=np.concatenate(histt5, axis=0) 
n4=np.concatenate(histt6, axis=0) 
n5=np.concatenate(histt7, axis=0)  
n6=np.concatenate(histt8, axis=0)

###PREPROCESSING pH BASED TICA DATA 
j=0
tica=[]                                                                                                                                                
for k in range(6):       
    f1=[]                                                                                                
    for i in range(int(len(tica_trajs)/3)):                                                                                              
        f=list(tica_trajs[3*i+2][j:j+800]);
        f1.append(f)                                                   
        print(3*i+2)                           
        j=j                                                                                                                                             
    f1=np.array(f1)                                                                                                                     
    da = np.concatenate(f1, axis=0)                                                                                 
    loop = range(len(da))                                                                                           
    str_list=list()                                                                                                                                                                                  
  
    da_int=da.astype(int)                                                                                               
    for x in loop:                                                                                                                                                                                   
                                                                                                                     
        a=str(da_int[x,:])                                                                                                                                            
        str_list.append(a)                                                                                                                                             
    unique=list(set(str_list))                                                                                         
    index=list()                                                                                                       
    number_unique = len(unique)                                                     
    print ("unique strings:  ", number_unique, "of total strings:", len(str_list))                                                                                                                   
                               
    tic=[]                                                                                                             
    for i in range(len(da)):                                                                                            
        tic.append(da[i])                       
    tica.append(tic)
    j=j+800

###CONCATENATION TICA DATA
tica=np.array(tica)
tica=np.concatenate(tica, axis=0)
print(np.shape(tica))

###QUANTIFICATION OF CONFORMATIONS
histtt2, newedgesX, newedgesY = np.histogram2d(tica[:,0], tica[:,1], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histtt3, newedgesX3, newedgesY3 = np.histogram2d(tica[:,0][0:8000], tica[:,1][0:8000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histtt4, newedgesX4, newedgesY4 = np.histogram2d(tica[:,0][8001:16000], tica[:,1][8001:16000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histtt5, newedgesX5, newedgesY5 = np.histogram2d(tica[:,0][16001:24000], tica[:,1][16001:24000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histtt6, newedgesX6, newedgesY6 = np.histogram2d(tica[:,0][24001:32000], tica[:,1][24001:32000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histtt7, newedgesX7, newedgesY7 = np.histogram2d(tica[:,0][32001:40000], tica[:,1][32001:40000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
histtt8, newedgesX8, newedgesY8 = np.histogram2d(tica[:,0][40001:48000], tica[:,1][40001:48000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])

o1=np.concatenate(histtt3, axis=0)
o2=np.concatenate(histtt4, axis=0) 
o3=np.concatenate(histtt5, axis=0) 
o4=np.concatenate(histtt6, axis=0) 
o5=np.concatenate(histtt7, axis=0)  
o6=np.concatenate(histtt8, axis=0)

###PREPROCESSING PH BASED TICA DATA 
j=0
tica=[]                                                                                                                                                
for k in range(6):                                                                                                                                                                                   
                                                                                                                      
    f1=[]                                                                                                
    for i in range(int(len(tica_trajs)/3)):                                                                                              
        f=list(tica_trajs[3*i+1][j:j+800])                      
        f1.append(f)                                                   
        print(3*i+1)                           
        j=j                                                                                                                                             
    f1=np.array(f1)                                                                                                                     
    da = np.concatenate(f1, axis=0)                                                                                 
    loop = range(len(da))                                                                                           
    str_list=list()                                                                                                                                                                                  
  
    da_int=da.astype(int)                                                                                               
    for x in loop:                                                                                                                                                                                   
                                                                                                                     
        a=str(da_int[x,:])                                                                                                                                            
        str_list.append(a)                                                                                                                                             
    unique=list(set(str_list))                                                                                         
    index=list()                                                                                                       
    number_unique = len(unique)                                                     
    print ("unique strings:  ", number_unique, "of total strings:", len(str_list))                                                                                                                   
                               
    tic=[]                                                                                                             
    for i in range(len(da)):                                                                                            
        tic.append(da[i])                       
    tica.append(tic)
    j=j+800

tica=np.array(tica)
tica=np.concatenate(tica, axis=0)
print(np.shape(tica))
###QUANTIFICATION OF CONFORMATIONS
hist2, newedgX, newedgY = np.histogram2d(tica[:,0], tica[:,1], bins=10)
hist3, newedgX3, newedgY3 = np.histogram2d(tica[:,0][0:8000], tica[:,1][0:8000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
hist4, newedgX4, newedgY4 = np.histogram2d(tica[:,0][8001:16000], tica[:,1][8001:16000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
hist5, newedgX5, newedgY5 = np.histogram2d(tica[:,0][16001:24000], tica[:,1][16001:24000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
hist6, newedgX6, newedgY6 = np.histogram2d(tica[:,0][24001:32000], tica[:,1][24001:32000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
hist7, newedgX7, newedgY7 = np.histogram2d(tica[:,0][32001:40000], tica[:,1][32001:40000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])
hist8, newedgX8, newedgY8 = np.histogram2d(tica[:,0][40001:48000], tica[:,1][40001:48000], bins=10, range=[[-4.9952, 0.935086], [-4.9281, 4.5184]])

m1=np.concatenate(hist3, axis=0)
m2=np.concatenate(hist4, axis=0) 
m3=np.concatenate(hist5, axis=0) 
m4=np.concatenate(hist6, axis=0) 
m5=np.concatenate(hist7, axis=0)  
m6=np.concatenate(hist8, axis=0)

###PLOTTING FUNCTION FOR FREE ENERGY FIGURE
fig=plt.figure()                           
fig.subplots_adjust(hspace=0.7, wspace=0.7) 
fig.text(0.5, 0.04, 'Temperature', ha='center')
fig.text(0.08, 0.5, 'Population', va='center', rotation='vertical')   
tr=[]                                      
for i in range(100):                                                                                     
    tr=[n1[i], n2[i], n3[i], n4[i], n5[i], n6[i]]
    row_based_plot_idx = row_based_idx(10, 10, i+1)
    ax=fig.add_subplot(10, 10, row_based_plot_idx)
    ax.plot(tr, color="blue", linewidth=2)  
    ax.xaxis.set_ticks(np.arange(0,6,1))
    matplotlib.rc('xtick', labelsize=7)
    matplotlib.rc('ytick', labelsize=7)
t=[]                                       
for i in range(100):                                                                                     
    t=[m1[i], m2[i], m3[i], m4[i], m5[i], m6[i]] 
    row_based_plot_idx = row_based_idx(10, 10, i+1)
    ax=fig.add_subplot(10, 10, row_based_plot_idx)
    ax.plot(t, color="red", linewidth=2)    
    ax.xaxis.set_ticks(np.arange(0,6,1))
    matplotlib.rc('xtick', labelsize=7)
    matplotlib.rc('ytick', labelsize=7)
r=[]                                       
for i in range(100):                                                                                     
    r=[o1[i], o2[i], o3[i], o4[i], o5[i], o6[i]] 
    row_based_plot_idx = row_based_idx(10, 10, i+1)
    ax=fig.add_subplot(10, 10, row_based_plot_idx)
    ax.plot(r, color="green", linewidth=2)  
    ax.xaxis.set_ticks(np.arange(0,6,1))
    matplotlib.rc('xtick', labelsize=7)
    matplotlib.rc('ytick', labelsize=7)
labels=["pH3", "pH5", "pH7"]
blue_patch = mpatches.Patch(facecolor='b', edgecolor='#000000')
red_patch = mpatches.Patch(facecolor='r', edgecolor='#000000')                                           
green_patch = mpatches.Patch(facecolor='g', edgecolor='#000000')
fig.legend(handles = [blue_patch, red_patch, green_patch], labels=labels, loc="upper right",borderaxespad=0.1)                                              
plt.subplots_adjust(right=0.95)
plt.show()

#TRACEPATH OF pH DEPENDENT CONFORMATIONS
import matplotlib.pyplot as pp
from seaborn.distributions import (_scipy_univariate_kde, _scipy_bivariate_kde)
cs=ax.contourf(X, Y, Z - Z.min(), cmap=plt.get_cmap("hsv"),levels=np.linspace(vmin, vmax, 10),alpha=1, vmin=vmin, vmax=vmax)

###EXTRACT FRAME FOR PH3
j=0
tica1=[]                                                                                                                                                
for k in range(6):       
    f1=[]                                                                                                
    for i in range(int(len(tica_trajs)/3)):                                                                                              
        f=list(tica_trajs[3*i][j:j+800])                      
        f1.append(f)                                                   
        print(3*i)                           
        j=j                                                                                                                                             
    f1=np.array(f1)                                                                                                                     
    da = np.concatenate(f1, axis=0)                                                                                 
    loop = range(len(da))                                                                                           
    str_list=list()                                                                                                                                                                                 
   
    da_int=da.astype(int)                                                                                               
    for x in loop:                                                                                                                                                                                  
                                                                                                                      
        a=str(da_int[x,:])                                                                                                                                            
        str_list.append(a)                                                                                                                                             
    unique=list(set(str_list))                                                                                         
    index=list()                                                                                                       
    number_unique = len(unique)                                                     
    print ("unique strings:  ", number_unique, "of total strings:", len(str_list))                                                                                                                  
                                
    tic=[]                                                                                                             
    for i in range(len(da)):                                                                                            
        tic.append(da[i])                       
    tica1.append(tic)
    j=j+800

###CONCATENATED ALL THE DATA INTO A SINGLE LIST 
tica1=np.array(tica1)
tica1=np.concatenate(tica1, axis=0)
print(np.shape(tica1))

###EXTRACT FRAME FOR PH7
j=0
tica2=[]                                                                                                                                                
for k in range(6):       
    f1=[]                                                                                                
    for i in range(int(len(tica_trajs)/3)):                                                                                              
        f=list(tica_trajs[3*i+2][j:j+800])                      
        f1.append(f)                                                   
        print(3*i+2)                           
        j=j                                                                                                                                             
    f1=np.array(f1)                                                                                                                     
    da = np.concatenate(f1, axis=0)                                                                                 
    loop = range(len(da))                                                                                           
    str_list=list()                                                                                                                                                                                 
   
    da_int=da.astype(int)                                                                                               
    for x in loop:                                                                                                                                                                                  
                                                                                                                      
        a=str(da_int[x,:])                                                                                                                                            
        str_list.append(a)                                                                                                                                             
    unique=list(set(str_list))                                                                                         
    index=list()                                                                                                       
    number_unique = len(unique)                                                     
    print ("unique strings:  ", number_unique, "of total strings:", len(str_list))                                                                                                                  
                                
    tic=[]                                                                                                             
    for i in range(len(da)):                                                                                            
        tic.append(da[i])                       
    tica2.append(tic)
    j=j+800

###CONCATENATED ALL THE DATA INTO A SINGLE LIST
tica2=np.array(tica2)
tica2=np.concatenate(tica2, axis=0)
print(np.shape(tica2))

###EXTRACT FRAME FOR PH5
j=0
tica3=[]                                                                                                                                                
for k in range(6):       
    f1=[]                                                                                                
    for i in range(int(len(tica_trajs)/3)):                                                                                              
        f=list(tica_trajs[3*i+1][j:j+800])                      
        f1.append(f)                                                   
        print(3*i+1)                           
        j=j                                                                                                                                             
    f1=np.array(f1)                                                                                                                     
    da = np.concatenate(f1, axis=0)                                                                                 
    loop = range(len(da))                                                                                           
    str_list=list()                                                                                                                                                                                 
   
    da_int=da.astype(int)                                                                                               
    for x in loop:                                                                                                                                                                                  
                                                                                                                      
        a=str(da_int[x,:])                                                                                                                                            
        str_list.append(a)                                                                                                                                             
    unique=list(set(str_list))                                                                                         
    index=list()                                                                                                       
    number_unique = len(unique)                                                     
    print ("unique strings:  ", number_unique, "of total strings:", len(str_list))                                                                                                                  
                                
    tic=[]                                                                                                             
    for i in range(len(da)):                      if __name__ == '__main__':
    main()                                                                      
        tic.append(da[i])                       
    tica3.append(tic)
    j=j+800

###CONCATENATED ALL THE DATA INTO A SINGLE LIST
tica3=np.array(tica3)
tica3=np.concatenate(tica3, axis=0)
print(np.shape(tica3))

###########MAIN FUNCTION FOR PLOT_TRACE2D
def plot_trace2d(data, obs=(0, 1), ts=1.0, cbar=True, ax=None, xlabel=None, ylabel=None, color=None, labelsize=14,cbar_kwargs=None, scatter_kwargs=None, plot_kwargs=None):
    if ax is None:
        ax = plt.gca()
    if scatter_kwargs is None:
        scatter_kwargs = {}
    if plot_kwargs is None:
        plot_kwargs = {}
    if not isinstance(obs, tuple):
        raise ValueError('obs must be a tuple')

    if isinstance(data, list):
        for item in data:
            prune = item[:, obs]
            ax.plot(prune[:, 0], prune[:, 1], **plot_kwargs)
    else:
        prune = data[:, obs]
        c = ax.scatter(prune[:, 0], prune[:, 1], c=color, **scatter_kwargs)

    if xlabel:
        ax.set_xlabel(xlabel, size=labelsize)
    if ylabel:
        ax.set_ylabel(ylabel, size=labelsize)
    return ax

##PLOTTING FUNCTION PLOT_TRACE2D FOR PH3
###TICKS THE X-AXIS AND Y-AXIS
xticks=np.arange(0.95, -4.5, -0.58)
yticks=np.arange(4.2, -4.9, -0.93)
ax=plt.axes()
cf=ax.contourf(X, Y, Z - Z.min(), cmap=plt.get_cmap("hsv"),levels=np.linspace(vmin, vmax, 10),alpha=0.5, vmin=vmin, vmax=vmax) 
ax.set(xlim=(-4.5, 0.95), ylim=(-4.9, 4.2), xlabel='tIC1', ylabel='tIC2')  
ax.set_xticks(xticks) 
ax.set_yticks(yticks)
ax.grid(True)
ax.grid(color='k', linestyle='-')  

###LABELLING OF THE CBAR
if cbar:
    if not cbar_kwargs:
        cbar_kwargs = {}
    if shade:
        mappable=cf
    else:
        mappable=cs
    plt.colorbar(mappable, **cbar_kwargs)
    
plot_trace2d(tica1, ts=1, ax=ax, scatter_kwargs={'s' :2}, color= "blue", xlabel='tIC1', ylabel='tIC2')
plt.show()

##PLOTTING FUNCTION PLOT_TRACE2D FOR PH5
###TICKS THE X-AXIS AND Y-AXIS
xticks=np.arange(0.95, -4.5, -0.58)
yticks=np.arange(4.2, -4.9, -0.93)
ax=plt.axes()
cf=ax.contourf(X, Y, Z - Z.min(), cmap=plt.get_cmap("hsv"),levels=np.linspace(vmin, vmax, 10),alpha=0.5, vmin=vmin, vmax=vmax) 
ax.set(xlim=(-4.5, 0.95), ylim=(-4.9, 4.2), xlabel='tIC1', ylabel='tIC2')  
###TICKS THE X-AXIS AND Y-AXIS
ax.set_xticks(xticks) 
ax.set_yticks(yticks)
ax.grid(True)
ax.grid(color='k', linestyle='-')  

###LABELLING OF THE CBAR
if cbar:
    if not cbar_kwargs:
        cbar_kwargs = {}
    if shade:
        mappable=cf
    else:
        mappable=cs
    plt.colorbar(mappable, **cbar_kwargs)
    
plot_trace2d(tica2, ts=1, ax=ax, scatter_kwargs={'s' :2}, color= "green", xlabel='tIC1', ylabel='tIC2')
plt.show()

##PLOTTING FUNCTION PLOT_TRACE2D FOR PH7
###TICKS THE X-AXIS AND Y-AXIS
xticks=np.arange(0.95, -4.5, -0.58)
yticks=np.arange(4.2, -4.9, -0.93)
ax=plt.axes()
cf=ax.contourf(X, Y, Z - Z.min(), cmap=plt.get_cmap("hsv"),levels=np.linspace(vmin, vmax, 10),alpha=0.5, vmin=vmin, vmax=vmax) 
ax.set(xlim=(-4.5, 0.95), ylim=(-4.9, 4.2), xlabel='tIC1', ylabel='tIC2')  
ax.set_xticks(xticks) 
ax.set_yticks(yticks)
ax.grid(True)
ax.grid(color='k', linestyle='-')  

###LABELLING OF THE CBAR
if cbar:
    if not cbar_kwargs:
        cbar_kwargs = {}
    if shade:
        mappable=cf
    else:
        mappable=cs
    plt.colorbar(mappable, **cbar_kwargs)
    
plot_trace2d(tica3, ts=1, ax=ax, scatter_kwargs={'s' :2}, color= "red", xlabel='tIC1', ylabel='tIC2')
plt.show()

if __name__ == '__main__':
    main()
