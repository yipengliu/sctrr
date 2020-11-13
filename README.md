# SCTRR

Code for our paper "Smooth Compact Tensor Ring Regression" in IEEE Transactions on Knowledge and Data Engineering.

Requirements
The algorithms have been implemented in MATLAB and based on Tensor Toolbox 2.6 (http://www.sandia.gov/~tgkolda/TensorToolbox/).


Contents
demo  --- a simulation test for the proposed SCTRR method
SCTRR --- the main function for the proposed SCTRR method
init_trmodel --- generate a random tensor rings with given size and rank
update-ul --- the solution for each subproblem with respect to U^(1), ..., U^(L)
updat_vm --- the solution for each subproblem with respect to U^(L+1), ..., U^(L+M)


Reference
 J. Liu, C. Zhu and Y. Liu, "Smooth Compact Tensor Ring Regression," in IEEE Transactions on Knowledge and Data Engineering, doi: 10.1109/TKDE.2020.3037131.
