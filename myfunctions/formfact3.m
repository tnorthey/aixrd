function[fq,fq3]=formfact3(an,q)
% Obtain form-factor from:
% http://lamp.tu-graz.ac.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

% input atom number (an) and vector q in Angstroms
% outputs form factor for atom number and q-range

iam_constants;  % H to F constants, if needed obtain more from website.

Nq=length(q);
fq=zeros(Nq,1);
for j=1:Nq
    for i=1:4
        fq(j)=fq(j)+aa(an,i)*exp(-bb(an,i)*(.25*q(j)/pi)^2);
    end
end
fq=fq+cc(an);

% Spherical 3D form-factor (could be useful later),
fq3=zeros(Nq,Nq,Nq);
for k1=1:Nq
    for k2=1:Nq
        for k3=1:Nq
            fq3(k1,k2,k3)=fq(k1);  % invariant with theta and phi!
        end
    end
end

% plot(q,fq);%

return