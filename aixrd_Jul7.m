% Author: Thomas Northey % Date: 11/6/14
% Description: Based on xftgen2_new.m

% Calculates analytical diffraction pattern for input file 'mldfile'
% with wavevector k0=(0,0,kz).
% Photon wavelength wl in angstroms. 
% Disregards GTO products with coefficients Mi*Mj*ci*cj <= cutoff.
% Outputs FT(theta,phi) with prec theta, phi increments.
% 7 July 2014: Current fastest version, with full analytical accuracy.

function[FT]=aixrd_Jul7(mldfile,prec,cutoff)
tic
prec=2*ceil(prec/2)+1;  % make prec odd, so FT(ph=0) is directly across from FT(ph=pi/2) on circle
% Read coefficients from Molden file
[bi,M,ga,c,l,m,n,xx,yy,zz,ppmo,moocc,~,Atoms]=mldread(mldfile);
% Disregard Mi*Mj*ci*cj <= cutoff
[I,J] = find_significant_gtos2(M,c,ppmo,moocc,cutoff);
bi=bi(I);Mi=M(I);Mj=M(J);ci=c(I);cj=c(J);li=l(I);lj=l(J);mi=m(I);mj=m(J);
ni=n(I);nj=n(J);xi=xx(I);xj=xx(J);yi=yy(I);yj=yy(J);zi=zz(I);zj=zz(J);
gi=ga(I);gj=ga(J);

Time =  num2str(ceil(length(I)*prec^2*5.745e-05));
disp(strcat('This will take ~',Time,' seconds (single 3.1 GHz core)'))
disp('(Lower precision or increase cutoff to decrease time.)')

% other definitions:
au2ang = 0.52917721092d0;
wl=1.4586382;                  % angstrom (8.5 keV)
lam=wl/au2ang;                 % wavelength in Angstroms
kk=2.d0*pi/lam;                % magnitude of k-vector 
nth=prec;nph=prec;
th=linspace(0.d0,1*pi,nth);    % 0 < theta < pi
ph=linspace(0.d0,2*pi,nph);    % 0 < phi < 2pi 
% keyboard
% ---------------------------------------------------------------------------------------
% analytical diffraction calculation (for MO range mo1-mo2):
disp('Transforming...')
FT=zeros(nth,nph);
for B = 1:nph-1                 % s(ph=0)=s(ph=2pi), no need to calc. twice
    for A = 1:nth        
        % s-vector definition (x-ray along z-axis):
        sx=kk*sin(th(A))*cos(ph(B)); sy=kk*sin(th(A))*sin(ph(B)); sz=kk*cos(th(A))-kk;
        for i = 1:length(I);                
            pref=bi(i)*Mi(i)*Mj(i)*ci(i)*cj(i);
            FT(A,B) = FT(A,B) + pref*ftswch3(gi(i),gj(i),li(i),lj(i),xi(i),xj(i),sx)*...
                                     ftswch3(gi(i),gj(i),mi(i),mj(i),yi(i),yj(i),sy)*...
                                     ftswch3(gi(i),gj(i),ni(i),nj(i),zi(i),zj(i),sz);            
        end
    end 
end
FT(:,nph)=FT(:,1);                  % FT(ph=0) = FT(ph=2pi)
Nelec0=sum(Atoms(:,2));
Nelec=max(max(abs(FT)));
disp(strcat('True Nelec=',num2str(Nelec0)));
disp(strcat('Integrated Nelec=',num2str(Nelec))); % to check
disp('Note: Integrated Nelec = True Nelec if cutoff=0');
% ---------------------------------------------------------------------------------------
% PLOT
figure('color','white');
polarplot3d(abs(FT).^2,'PolarGrid',{1 1},'plottype','contour','ContourLines',200,...
    'TickSpacing',45,'MeshScale',[.1 .1]);
colorbar;
toc
% ---------------------------------------------------------------------------------------

return
