%% _Ab initio_ molecular diffraction.
% 
% T. Northey, August 2015.
%
% School of Chemistry, University of Edinburgh, West Mains Road, 
% Edinburgh EH9 3FJ, United Kingdom.
%
% Reference:
% T. Northey, N. Zotev, A. Kirrander, J. Chem. Theory Comput. 
% (2014) DOI: 10.1021/ct500096r
%
%% Description
% Reads molden file (generated with molpro) from input folder
% and calculates the molecular form-factor for elastic x-ray diffraction,
% using settings specified in setup2.txt (provided),
% 
%  Example setup2.txt:
% 
%  ==================================
%  title  Nq  wl  dim  method
% 
%  h2 49 1.5 sph ai
%  ethylene 99 1.5 x iam
% 
%  ==================================
% 
% Reads from input/h2.mld and calculates for momentum transfer vector (q) 
% of length 49, wavelength 1.5 bohr, where q=(4*pi/wl)*sin(theta/2), 
% and theta is the scattering vector theta=[0;pi], 'sph' dimensions option 
% (spherical coordinates), for 3D output f(|q|,theta,phi), where 
% size(f)=[49,49,49], other options are 'x','y','z' for 1D outputs f(qx), 
% f(qy), or f(qz), and 'ai' option (ab initio) for ab initio x-ray 
% diffraction method, the other option is 'iam' method (Independent Atom 
% Model). Then, results are saved in results folder as .mat file.
% 
% As there is a second line in the example setup2.txt, the code will then 
% perform a second calculation based on input/ethylene.mld with Nq=49, 
% wl=1.5, dim='x' and the IAM method, size(f)=[99,1].
% 
%% Usage
% Modify setup2.txt with chosen options, make sure the correct molden
% file is in input folder. Then run this function.
%
%% Efficiency
% The efficiency for the AIXRD method scales with the square of the total 
% number of Gaussians that make up the wave function. Large molecules and 
% basis sets can take a very long time. There is a time estimator within 
% the function. 
%
%% Theory
% The theory is shown in the reference paper.
%% Function begins
function [] = AIXRD2016()

% Calculates elastic x-ray diffraction from a list file, 'setup.txt'.
addpath(strcat(pwd,'/myfunctions'))
fid = fopen('setup2.txt', 'r'); % Open file for reading                                              
tline = fgetl(fid); % Read line from file, removing newline characters
tline = fgetl(fid); % Read line from file, removing newline characters
c=0;
while ~feof(fid)
    c=c+1;
    tline = fgetl(fid); % Read line from file, removing newline characters
    tmp=textscan(tline,'%s %f %f %s %s');
    mldfile(c)=strcat('input/',tmp{1},'.mld');
    Nq(c)=tmp{2};
    wl(c)=tmp{3};
    dim(c)=tmp{4};
    method(c)=tmp{5};
end

for i=1:c
    aixrd_Aug2015(method{i},dim{i},mldfile{i},Nq(i),wl(i));
end

return

%% Scattering calculation
function [] = aixrd_Aug2015(method,dim,mldfile,Nq,wl)
%=======================================
tic  % start timer
[~,rtitle,~] = fileparts(mldfile);    % define title from mldfile.
rtitle=strcat(rtitle,'_',method,'_',dim,'_Nq',num2str(Nq),'_wl',num2str(wl));

%========================================
% calculate Fq
if strcmpi(method,'iam')
    [Fq,Atoms,C]=AIXRD2016_calcFq_iam(mldfile,Nq,wl,dim);
elseif strcmpi(method,'ai')
    [Fq,Atoms,C]=AIXRD2016_calcFq_ai(mldfile,Nq,wl,dim);
end
Iq=abs(Fq).^2;     % Intensity is form-factor absolute squared ff*
%================================================================================
% ROTATIONAL AVERAGE (if using spherical coordinates in q)
if strcmp(dim,'sph')
    Irot=rotavg0(C{2},C{3},Iq);  % ensemble-average (each particle scatters incoherently)
    Ideb=iamcalc2(Atoms(:,2:end),C{1});  % Debye approx.
    Frot=rotavg0(C{2},C{3},Fq);  % particle scatters coherently (with itself usually)
end
%================================================================================
% SAVE to output
% save output variables in .mat file
matname=strcat('results/',rtitle,'.mat');
if strcmp(dim,'sph') 
    save(matname,'rtitle','Fq','C','Frot','Irot','Ideb');
elseif strfind(dim,'det')
    save(matname,'rtitle','Fq','C');
else
    save(matname,'rtitle','Fq','C');
end
%================================================================================
toc  % stop timer
disp('Done.')
return % function end
