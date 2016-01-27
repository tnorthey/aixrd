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
function [] = AIXRD2015_v3pt1_stable()

% Calculates elastic x-ray diffraction from a list file, 'setup.txt'.

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

Nq=2*ceil(Nq/2)-1;  % want Nq to be odd so q=0 is in linspace(-qmax,qmax,Nq)
                    % because q=0 is important (central maximum)

k0=2*pi/wl;         % incoming vector magnitude |k|
qmax=2*k0;          % q = k0 - k, where |k0|=|k| (elastic), so qmax=2|k|
if Nq==1
    q=0;
else
    q=linspace(0,qmax,Nq);  % define momentum transfer vector q
end

% coordinate choice (dim='x','y','z','xy','xz','yz','xyz','sph')
[qx,qy,qz,Fq,C] = choosecoord2(dim,q);

if size(C,1)==3
    r=C{1};
    th=C{2};
    ph=C{3};
end

Q=cell(3,1);  % make convenient cell arrays
Qe=Q; Qi=Q;
Q{1}=qx; Q{2}=qy; Q{3}=qz;
for k=1:3
    Qe{k}=-.25*Q{k}.^2;
    Qi{k}=.5i*Q{k};
end

%========================================
% calculate Fq
[Fat,Fmol,Atoms]=AIXRD2015_calcFq_v3pt1(method,dim,mldfile,Nq,Q,Qe,Qi,Fq,q);
Fq=Fat+Fmol;
Iq=abs(Fq).^2;     % Intensity is form-factor absolute squared ff*
%================================================================================
% ROTATIONAL AVERAGE (if using spherical coordinates in q)
if strcmp(dim,'sph')
    Irot=rotavg0(th,ph,Iq);  % ensemble-average (each particle scatters incoherently)
    Ideb=iamcalc2(Atoms(:,2:end),q);  % Debye approx.
    Frot=rotavg0(th,ph,Fq);  % particle scatters coherently (with itself usually)
end
%================================================================================
% SAVE to output
% save output variables in .mat file
matname=strcat('results/',rtitle,'.mat');
if strcmp(dim,'sph') 
    save(matname,'rtitle','Fat','Fmol','r','th','ph','Frot','Irot','Ideb');
elseif strfind(dim,'det')
    save(matname,'rtitle','Fat','Fmol','r','th','ph');
else
    save(matname,'rtitle','Fat','Fmol','q');
end
%================================================================================
toc  % stop timer
disp('Done.')
return % function end
