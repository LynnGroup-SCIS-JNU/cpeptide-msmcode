# **C-Peptide-MSMCode**  
## Temperature and pH dependent MD simulation Analysis of C-peptide.

![Alt text](Results/C-peptide.jpg?raw=true "Title")

### To predict pH-dependent transition of conformations at various temperature using Markov state models on a standard free energy surface 
The code splits into 5 sections 
1. Data preprocessing 
2. Fel construction 
3. Population quantification 
4. pH dependent conformation distribution 
5. Tracing path of trajectories 

Execute the cpeptide.py to perform the analysis and find out the distribution of congormation on FEL.
```
python cpeptide.py
```

## *Package content* :
- Dataset :  trajectories *.nc & intial conformation s.pdb.
- Results : msm_mdl.pkl & *.png
- Analysis code : cpeptide.py

## *Dependencies* :
- [MSMBuilder v3.8](http://msmbuilder.org/3.8.0/)
- [MdTraj v1.9.4](https://mdtraj.org/1.9.4/index.html)
- [conda](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)
- python v3.7 
- python libraries

## *Installation* :
```
 bash Miniconda3-latest-Linux-x86_64.sh
 conda create -n myenv python=3.7
 source activate myenv
 conda install -c omnia msmbuilder
 conda install -c conda-forge mdtraj
 conda install package-name
```

## *Contact Information* :
For further details and bugs feel free to write  
- *Ruhar*,  ruhar63_sit@jnu.ac.in 
- *Andrew Lynn*, andrew@jnu.ac.in
