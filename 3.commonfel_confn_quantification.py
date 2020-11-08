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

