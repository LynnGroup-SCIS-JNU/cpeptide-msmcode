![Alt text](Results/C-peptide.jpg?raw=true "Title")

# **C-Peptide-MSMCode**  
## Temperature and pH dependent MD simulation Analysis of C-peptide.

pH dependent state transition of c-peptide and evaluation of conformations on FEL at various temperature.

The single code is splitted into 5 sections for detailed analysis which includes
1. Data preprocessing 
2. Common fel construction of all pH and temperatu:re simulated trajectories
3. Quantification of conformation on fel
4. Quantification of pH dependent conformations on FEL
5. Tracing path of pH dependent trajectories on FEL. 

Execute the cpeptide.py to perform the analysis and find out the distribution of congormation on FEL.
- python cpeptide.py

## Package content 
- Dataset : includes trajectories *.nc and intial conformation s.pdb.
- Results 
- Analysis code 

## *Dependencies* 
- [MSMBuilder v3.8] (http://msmbuilder.org/3.8.0/)
- MdTraj v1.9.4
- conda
- python v3.7 
- python libraries

## Installation
- bash Anaconda-latest-Linux-x86_64.sh
- conda install -c omnia msmbuilder
- conda install -c conda-forge mdtraj
- conda create -n myenv python=3.7
- source activate myenv
- conda install package-name

## *Reference for the previous information encouraged to current analysis*
- https://www.biorxiv.org/content/10.1101/084020v2.full 
- https://pubs.acs.org/doi/10.1021/jp401587e  
- https://github.com/msmbuilder/paper
- http://scaledmd.ucsd.edu/

## Contact Information
For further details and bugs feel free to write  
- *Ruhar*,  ruhar63_sit@jnu.ac.in 
- *Andrew Lynn*, andrew@jnu.ac.in
