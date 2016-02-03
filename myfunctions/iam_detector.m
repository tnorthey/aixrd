function [Fq] = iam_detector(Atoms,dim,q)
%IAM_DETECTOR Summary of this function goes here
%   Detailed explanation goes here

Nat=Atoms(end,1);   % no. of atoms
R0=Atoms(:,3:5);    % nuclear geometry (Bohr)
an=Atoms(:,2);      % atomic numbers
% au2ang = 0.52917721092d0; % convert Ang coordinates to au
Nq=length(q);
kv=max(q)/2;        % scattered vector k
th=linspace(0,pi,Nq);
ph=linspace(0,2*pi,Nq);

Fq=zeros(Nq,Nq);
%========================================
% calculate Fq
% Independent atom model method,
% total 3D form-factor, f(q),
for i = 1:Nat           % loop through atoms         
    for j=1:Nq          % theta
        for k=1:Nq      % phi
            if strcmp(dim,'detz')
                qx=kv*sin(th(j))*cos(ph(k));
                qy=kv*sin(th(j))*sin(ph(k));
                qz=kv*cos(th(j))-kv;
            elseif strcmp(dim,'dety')
                qz=kv*sin(th(j))*cos(ph(k));
                qx=kv*sin(th(j))*sin(ph(k));
                qy=kv*cos(th(j))-kv;
            elseif strcmp(dim,'detx')  
                qy=kv*sin(th(j))*cos(ph(k));
                qz=kv*sin(th(j))*sin(ph(k));
                qx=kv*cos(th(j))-kv;                
            end
            fq=formfact(an(i),sqrt(qx^2+qy^2+qz^2));            
            phs=exp(1i*(qx*R0(i,1)+qy*R0(i,2)+qz*R0(i,3)));
            Fq(j,k) = Fq(j,k) + fq*phs;    
        end
    end
end
%========================================
return

function[fq]=formfact(an,q)
% Obtain form-factor from:
% http://lamp.tu-graz.ac.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

% input atom number (an) and vector q in Angstroms
% outputs form factor for atom number and q-range

au2ang = 0.52917721092d0; % convert Ang coordinates to au at the end
iam_constants;  % H to F constants, if needed obtain more from website.
q=q/au2ang;

Nq=length(q);
fq=zeros(Nq,1);
for j=1:Nq
    for i=1:4
        fq(j)=fq(j)+aa(an,i)*exp(-bb(an,i)*(.25*q(j)/pi)^2);
    end
end
fq=fq+cc(an);

% plot(q,fq);%

return

