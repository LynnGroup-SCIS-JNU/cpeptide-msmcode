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

