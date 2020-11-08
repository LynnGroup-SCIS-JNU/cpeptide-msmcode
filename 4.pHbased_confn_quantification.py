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

