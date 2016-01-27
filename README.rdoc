How to run AIXRD
T. Northey 2015
Main code: AIXRD_v3_stable.m

Usage:
Modify setup2.txt, containing 5 columns, title, Nq, wl, dim, method.
title: Arbitrary title, output .mat file will contain this title
Nq: Number of points in q-vector 
wl: x-ray wavelength in atomic units (Bohr)
dim: dimension, can be x,y,z which gives 1D slices through q of size (Nq,1);
detx, dety, detz, which give 2D detector images (functions of scattering 
angles theta, phi) of size (Nq,Nq) for incoming x-rays from the x, y, and z 
direction respectively; xyz, which gives a 3D grid in q of size (Nq,Nq,Nq); 
finally sph, which gives a 3D grid in q in spherical coordinates of size 
(Nq,Nq,Nq), and by default gives rotational-average and Debye approximation.
method: Either iam or ai, which means independent atom model (IAM, v. fast) 
approximation, or ab initio x-ray diffraction (AIXRD, slow) from [1]. 

setup2.txt looks like:
%-------------------------------%
title  Nq  wl  dim  method

hehe 49 3 detz ai
%-------------------------------%

In this case, running AIXRD_v3_stable.m will use data from the file 
input/hehe.mld (dihelium molden file from molpro).
Appropriate results get saved to results/ folder in .mat format.

Reference:
[1]. Thomas Northey, Nikola Zotev, and Adam Kirrander
School of Chemistry, University of Edinburgh, West Mains Road, Edinburgh EH9 3JJ, United Kingdom
J. Chem. Theory Comput., 2014, 10 (11), pp 4911–4920
DOI: 10.1021/ct500096r
Publication Date (Web): October 10, 2014