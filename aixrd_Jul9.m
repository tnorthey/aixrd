function [Fs] = aixrd_Jul9(mldfile,wl,Ns,cutoff)
% T. Northey 9/7/14
% Function calculates f(s) = FT[rho_ab-initio(r)](s) on an Ns^3 grid 
% for -smax <= s <= smax where smax=4*pi/wl, 
% rho_ab-initio(r) is general ab-intio electron density. 
% Units for wavelength, wl not specified. 
% Must select same units as molecular geometry units (au or angstrom)!
% Uses Molpro output file in Molden format, 'mldfile'.

% Additional notes:

% Improved from aixrd_Jul7 as instead of calculating f(theta,phi) for 
% a particular x-ray direction, the much more general f(s) is calculated. 
% Afterwards incoming x-ray direction can be considered.

% Takes advantage of f(s) = f(sx)*f(sy)*f(sz) by only calculating 
% along each axis and then multiplying. This is different to aixrd_Jul7.

% As such, calculation time scales as ~ Ngto^2 * Ns^3
% wheras aixrd_Jul7 scales as ~ 3*Ngto^2 * Ns^2.
% In general, this code is faster for Ns < 150.
% (It is faster because there are less calls to ftswch2.)

%=======================================

tic
[bi,M,ga,c,l,m,n,xx,yy,zz,ppmo,moocc,~,Atoms] = mldread(mldfile);           % Read necessary coefficients from molden file
[I,J] = find_significant_gtos2(M,c,ppmo,moocc,cutoff);                      % Find indices of all GTO products with pre-coeffs > cutoff 
bi=bi(I);Mi=M(I);Mj=M(J);ci=c(I);cj=c(J);li=l(I);lj=l(J);mi=m(I);mj=m(J);
ni=n(I);nj=n(J);xi=xx(I);xj=xx(J);yi=yy(I);yj=yy(J);zi=zz(I);zj=zz(J);
gi=ga(I);gj=ga(J);

Ns=2*ceil(Ns/2)-1;  % want Ns to be odd so s=0 is in linspace(-smax,smax,Ns)    
                    % because s=0 is important (centre of detector)

Time =  num2str(ceil(3*length(I)*Ns*2.4084e-05+Ns^3*length(I)*2.8334e-07)); % An estimate for the time a calculation will take.
disp(strcat('This will take ~',Time,' seconds (single 3.1 GHz core)'))
disp('(Lower precision to decrease time.)')

Fsx=zeros(Ns,length(I)); Fsy=Fsx; Fsz=Fsx;                                  % f(sx), f(sy), f(sz)

k0=2*pi/wl;                                                                 % incoming vector magnitude |k|
smax=2*k0;                                                                  % s = k0 - k, where |k0|=|k| (elastic), so smax=2|k|
ss=linspace(-smax,smax,Ns);                                                 % define momentum transfer vector s
Fs=zeros(Ns,Ns,Ns);
for k=1:Ns               % loop through s
    for i = 1:length(I)          % loop through significant GTOs  
 
        % F=F(sk,gij), FT of Gaussian product Fij at sk is F(k,i)
        pref=(bi(i)*Mi(i)*Mj(i)*ci(i)*cj(i))^(1/3);
        Fsx(k,i) = pref*ftswch2(gi(i),gj(i),li(i),lj(i),xi(i),xj(i),ss(k));
        Fsy(k,i) = pref*ftswch2(gi(i),gj(i),mi(i),mj(i),yi(i),yj(i),ss(k));
        Fsz(k,i) = pref*ftswch2(gi(i),gj(i),ni(i),nj(i),zi(i),zj(i),ss(k));
    end 
end
%========================================
%========================================
% complete Fs
for i = 1:length(I)         % loop through significant GTOs
    for k1=1:Ns
        for k2=1:Ns       
            for k3=1:Ns       
                Fs(k1,k2,k3) = Fs(k1,k2,k3) + Fsx(k1,i)*Fsy(k2,i)*Fsz(k3,i);
            end
        end
    end
end
Is=abs(Fs).^2;              % Intensity is form-factor absolute squared ff*
%========================================
      
% For checking purposes,
maxFs=max(max(max(abs(Fs)))); % f(s=0) will equal Nelec
Nelec=maxFs;
Nelec0=sum(Atoms(:,2));       % Nelec from molden file
disp(strcat('True Nelec=',num2str(Nelec0)));
disp(strcat('Integrated Nelec=',num2str(Nelec))); % to check
% disp('Note: Integrated Nelec = True Nelec if cutoff=0');

% PLOT
fig1=slice(ss,ss,ss,Is,0,0,0);  % cannot plot 3-D function easily, choose to plot intensity slices through sx=sy=sz=0
set(fig1,'edgecolor','none')
title('I(q)') 
xlabel('q_x (au^{-1})');ylabel('q_y (au^{-1})');zlabel('q_z (au^{-1})');
colorbar;

savefig('aixrd_output.fig');			% save slice figure

save('aixrd_output','mldfile','Is','Fs','ss');	% save output variables in .mat file

toc
disp('Done.')
% keyboard
% function ends
return

