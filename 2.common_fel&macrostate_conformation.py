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
