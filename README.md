# **cpeptide-msmcode**

#README for pH dependent state transition of c-peptide and evaluation of conformations on FEL at various temperature.
For further details and bugs contact Ruhar, PhD@JNU and Andrew Lynn, Prof@JNU. 
(ruhar63_sit@jnu.ac.in and andrew@jnu.ac.in)

This folder includes dataset, example results and analysis code for temperature and pH dependent MD simulation of C-peptide.
The trajectory are *.nc and initial conformation is s.pdb.

## *Supporting softwares* :
- MSMBuilder
- MdTraj
- conda
- python v3.7 
- python libraries

## *Reference for the previous information encouraged to current analysis*
- 1.https://www.biorxiv.org/content/10.1101/084020v2.full or https://github.com/msmbuilder/paper
- 2.https://pubs.acs.org/doi/10.1021/jp401587e or http://scaledmd.ucsd.edu/
 
Execute the cpeptide.py to perform the analysis and find out the distribution of congormation on FEL.

The single code is splitted into 5 sections for detailed analysis which includes
1. data preprocessing 
2. common fel construction of all pH and temperature simulated trajectories
3. quantification of conformation on fel
4. quantification of pH dependent conformations on FEL
5. tracing path of pH dependent trajectories on FEL. 