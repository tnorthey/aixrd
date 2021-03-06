# Ab Initio X-ray Diffraction

## Usage

Main code: AIXRD2016.m

Modify setup2.txt, which contains 5 columns, title, Nq, wl, dim, method. 
 
*title* (string): Arbitrary, results/'output'.mat will be created with 'output' based on title.  
*Nq* (integer): Number of points in momentum transfer vector (q, inverse bohr).  
*wl* (float): X-ray wavelength in atomic units (Bohr).  
*dim* (string): Dimension. Can either be x, y, or z which gives 1D slices through q of size (Nq,1); detx, dety, detz, which give 2D detector images (functions of scattering angles theta, phi) of size (Nq,Nq) for incoming x-rays from the x, y, and z direction respectively; xyz, which gives a 3D grid in q of size (Nq,Nq,Nq); finally sph, which gives a 3D grid in q in spherical coordinates of size (Nq,Nq,Nq), and by default gives rotationally-averaged curve and Debye approximation (for comparison).  
*method* (string): Either iam or ai. That is, independent atom model (IAM), or ab initio x-ray diffraction (AIXRD) [1]. 

#### setup2.txt format:

% beginning of file   
title  Nq  wl  dim  method

title1 Nq1 wl1 dim1 method1  
title2 Nq2 wl2 dim2 method2  
% end of file 

Run AIXRD2016.m to iterate over setup2.txt.  
Results get saved to results/ folder in .mat format.

## Reference

1. Thomas Northey, Nikola Zotev, and Adam Kirrander, J. Chem. Theory Comput., 2014, 10 (11), 4911–4920.  
DOI: [10.1021/ct500096r](http://dx.doi.org/10.1021/ct500096r)

---
