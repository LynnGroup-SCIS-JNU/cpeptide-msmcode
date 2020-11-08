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
