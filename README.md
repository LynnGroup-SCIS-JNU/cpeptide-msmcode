# **C-Peptide-MSMCode**  
## Temperature and pH dependent MD simulation Analysis of C-peptide.

![Alt text](Results/C-peptide.jpg?raw=true "Title")

pH dependent state transition of c-peptide and evaluation of conformations on FEL at various temperature.

The single code is splitted into 5 sections for detailed analysis which includes
1. Data preprocessing 
2. Fel construction for all pH and temperature simulated trajectories
3. Population quantification on fel
4. pH dependent distribution of conformation on FEL
5. Tracing path of pH dependent trajectories on FEL 

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
 conda install -c omnia msmbuilder
 conda install -c conda-forge mdtraj
 conda create -n myenv python=3.7
 source activate myenv
 conda install package-name
```

## *Contact Information* :
For further details and bugs feel free to write  
- *Ruhar*,  ruhar63_sit@jnu.ac.in 
- *Andrew Lynn*, andrew@jnu.ac.in
