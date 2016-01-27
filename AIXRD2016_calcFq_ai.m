function[Fq,Atoms,C]=AIXRD2016_calcFq_ai(mldfile,Nq,wl,dim)

% Setup q-array,
[C,Q,Qe,Qi,Fq,~] = setupq(Nq,wl,dim);

% Read Atom information from molden file,
Atoms = Atomsread(mldfile);

%========================================
% calculate Fq
Fii=zeros(size(Fq));
Fij=Fii;

% Ab initio x-ray diffraction method,
disp('Using ab initio method...')
[b,M,ga,c,l,m,n,xx,yy,zz,ppmo,moocc,~,~,~,~,~] = mldread_g(mldfile,0);    
cutoff=1e-9;    
Ng=length(c); % number of primitives

li=l;mi=m;ni=n;
% 3D GTO coordinates,
ri=[xx yy zz];

%========================================
% GTO normalisation,
% Normalisation factors in front of each prim. Gaussian
normA=(2/pi)^0.75;
normB=2.^(li+mi+ni);
normC=ga.^((2*li+2*mi+2*ni+3)/4);
normD=factd(2*li-ones(size(li))) .* factd(2*mi-ones(size(mi))) .* factd(2*ni-ones(size(ni)));
nrmi= normA .* normB .* normC ./ normD.^ 0.5;

gm = 2*ga;   % exponent addition
gu = gm.^-1;     

% pre-coeffs,
pre0=(pi*gu).^1.5.*b.*nrmi.^2.*M.^2.*c.^2;

for i=1:Ng     % diagonal terms i=j

    if pre0(i)>0   % only calculate for non-zero terms

        ex1=ones(size(Fq));                         % exponential term

        % 3D vector for convenience,
        Li=[li(i) mi(i) ni(i)];

        for j=1:3                                   % loop over x,y,z / l,m,n           

            ply=1;                                  % Li=Lj=0, polynomial term    

            if Li(j)==1
                A = Qi{j}*gu(i)+ri(i,j);
                Ri=ri(i,j);
                ply = A.^2+.5*gu(i)-2*Ri*A+Ri^2;               
            elseif Li(j)==2    
                A = Qi{j}*gu(i)+ri(i,j);
                Ri=ri(i,j);
                ply = A.^4-4*Ri*A.^3+...
                (3*gu(i)+2*Ri^2+8*Ri)*A.^2-...
                (2*Ri^3+2*Ri^3+6*Ri*gu(i))*A+...
                0.5*(Ri^2+Ri^2+4*Ri^2)*gu(i)+...
                0.75*gu(i)^2+Ri^4;               
            elseif Li(j)>0    
                ply=zeros(size(Fq));
                Ri=ri(i,j);
                for mm=0:2*Li(j)
                    v=2*Li(j)-mm;
                    for p=0:v
                        if mod(v-p,2) == 0 % if even   
                            ply = ply +...
                            nchoosek(2*Li(j),mm)*nchoosek(v,p)...
                            *(-Ri)^mm*(Ri+Qi{j}*gu(i)).^p...
                            *factd(v-p-1)*((.5*gu(i))^(.5*(v-p)));
                        end
                    end
                end                     
            end    % end if

            % exponential terms; Gaussian and phase,
            ex1=ex1.*ply.*exp(Qe{j}*gu(i)+1i*Q{j}*ri(i,j));

        end

        % a special case of atomic terms (i=j),
        Fii = Fii + pre0(i)*ex1;  
    end

end % end loop over Gaussian primitives

% when i not equal to j...
% Find indices of GTO products with pre-coeffs > cutoff,
[I,J] = find_significant_gtos3(M,c,ga,xx,yy,zz,ppmo,1,moocc,cutoff);
% make arrays correspinding to these indices,
bi=b(I);Mi=M(I);Mj=M(J);ci=c(I);cj=c(J);li=l(I);lj=l(J);mi=m(I);mj=m(J);
ni=n(I);nj=n(J);xi=xx(I);xj=xx(J);yi=yy(I);yj=yy(J);zi=zz(I);zj=zz(J);
gi=ga(I);gj=ga(J);

% 3D GTO coordinates,
ri=[xi yi zi];
rj=[xj yj zj];

Ngp=length(I);      % no. of GTO products to calculate

% Estimate approx. calculation time,                    
if strcmp(dim,'x') || strcmp(dim,'y') || strcmp(dim,'z')
    Time = num2str(round(3*Ngp*Nq*3e-08));
elseif strcmp(dim,'xy') || strcmp(dim,'xz') || strcmp(dim,'yz') ...
        || strcmp(dim,'detx') || strcmp(dim,'dety') || strcmp(dim,'detz')
    Time = num2str(round(3*Ngp*Nq^2*3e-08));    
elseif strcmp(dim,'xyz') || strcmp(dim,'sph')
    Time = num2str(round(3*Ngp*Nq^3*3e-08));        
end
disp(strcat('This will take ~',Time,' seconds (single 3.1 GHz core)'))
disp('(Lower precision to decrease time.)')

%========================================
% GTO normalisation,
% Normalisation factors in front of each prim. Gaussian
normA=(2/pi)^0.75;
normB=2.^(li+mi+ni);
normC=gi.^((2*li+2*mi+2*ni+3)/4);
normD=factd(2*li-ones(size(li))) .* factd(2*mi-ones(size(mi))) .* factd(2*ni-ones(size(ni)));
nrmi= normA .* normB .* normC ./ normD.^ 0.5;

normB=2.^(lj+mj+nj);
normC=gj.^((2*lj+2*mj+2*nj+3)/4);
normD=factd(2*lj-ones(size(lj))) .* factd(2*mj-ones(size(mj))) .* factd(2*nj-ones(size(nj)));
nrmj= normA .* normB .* normC ./ normD.^ 0.5;
%========================================

gm = gi+gj;   % exponent addition
gu = gm.^-1;     

% GTO overlap,
K=exp(-gi.*gj.*gu.*((xi-xj).^2+(yi-yj).^2+(zi-zj).^2));

% pre-coeffs,
pre0=(pi*gu).^1.5.*bi.*nrmi.*nrmj.*Mi.*Mj.*ci.*cj.*K;

rij=zeros(length(li),3);
for j=1:3
    rij(:,j) = (gi.*ri(:,j)+gj.*rj(:,j)).*gu;   % weighted GTO centre
end

for i = 1:Ngp                   % loop through significant GTO products    

    ex1=ones(size(Fq));                         % exponential term

    % 3D vector for convenience,
    Li=[li(i) mi(i) ni(i)];
    Lj=[lj(i) mj(i) nj(i)];

    for j=1:3                                   % loop over x,y,z / l,m,n           

        ply=1;                                  % Li=Lj=0, polynomial term    

        if Li(j)==0 && Lj(j)==1                 % higher order terms
            ply=Qi{j}*gu(i)+rij(i,j)-rj(i,j);
        elseif Li(j)==1 && Lj(j)==0    
            ply=Qi{j}*gu(i)+rij(i,j)-ri(i,j);            
        elseif Li(j)==1 && Lj(j)==1
            A = Qi{j}*gu(i)+rij(i,j);
            Ri=ri(i,j); Rj=rj(i,j);
            ply = A.^2+.5*gu(i)-(Ri+Rj)*A+Ri*Rj;               
        elseif Li(j)==0 && Lj(j)==2
            A = Qi{j}*gu(i)+rij(i,j);
            Rj=rj(i,j);
            ply = A.^2+0.5*gu(i)-2*Rj*A+Rj^2;  
        elseif Li(j)==2 && Lj(j)==0    
            A = Qi{j}*gu(i)+rij(i,j);
            Ri=ri(i,j);
            ply = A.^2+0.5*gu(i)-2*Ri*A+Ri^2;
        elseif Li(j)==1 && Lj(j)==2 
            A = Qi{j}*gu(i)+rij(i,j);
            Ri=ri(i,j); Rj=rj(i,j);
            ply = A.^3-(2*Rj+Ri)*A.^2+...
            (Rj^2+2*Rj*Ri+1.5*gu(i))*A-...
            (2*Rj+Ri)*0.5*gu(i)-Rj^2*Ri;
        elseif Li(j)==2 && Lj(j)==1    
            A = Qi{j}*gu(i)+rij(i,j);
            Ri=ri(i,j); Rj=rj(i,j);
            ply = A.^3-(2*Ri+Rj)*A.^2+...
            (Ri^2+2*Ri*Rj+1.5*gu(i))*A-...
            (2*Ri+Rj)*0.5*gu(i)-Ri^2*Rj;
        elseif Li(j)==2 && Lj(j)==2    
            A = Qi{j}*gu(i)+rij(i,j);
            Ri=ri(i,j); Rj=rj(i,j);
            ply = A.^4-2*(Ri+Rj)*A.^3+...
            (3*gu(i)+Ri^2+Rj^2+4*Ri*Rj)*A.^2-...
            (2*Ri^2*Rj+2*Ri*Rj^2+3*(Ri+Rj)*gu(i))*A+...
            0.5*(Ri^2+Rj^2+4*Ri*Rj)*gu(i)+...
            0.75*gu(i)^2+Ri^2*Rj^2;               
        elseif Li(j)>2 || Lj(j)>2    
            ply=zeros(size(Fq));
            Ri=ri(i,j); Rj=rj(i,j); Rij=rij(i,j);
            for m=0:Li(j)
                for n=0:Lj(j)
                    v=Li(j)+Lj(j)-(m+n);
                    for p=0:v
                        if mod(v-p,2) == 0 % if even   
                            ply = ply +...
                            nchoosek(Li(j),m)*nchoosek(Lj(j),n)*nchoosek(v,p)...
                            *(-Ri)^m*(-Rj)^n*(Rij+Qi{j}*gu(i)).^p...
                            *factd(v-p-1)*((.5*gu(i))^(.5*(v-p)));
                        end
                    end
                end
            end            
        end    % end if

        % exponential terms; Gaussian and phase,
        ex1=ex1.*ply.*exp(Qe{j}*gu(i)+1i*Q{j}*rij(i,j));

    end

    if ri(i,:)==rj(i,:)
        Fii = Fii + 2*pre0(i)*ex1;    % atomic terms
    else
        Fij = Fij + 2*pre0(i)*ex1;  % molecular terms
    end
end % end loop over Gaussian products

%========================================
% For checking purposes,
Fq=Fii+Fij;   % add diagonal and non-diagonal terms
maxFq=max(max(max(abs(Fq)))); % fmax = f(q=0) = Nelec
Nelec=maxFq;
Nelec0=sum(Atoms(:,2));       % Nelec from molden file
disp(strcat('True Nelec=',num2str(Nelec0)));
disp(strcat('Integrated Nelec=',num2str(Nelec))); % to check

return
