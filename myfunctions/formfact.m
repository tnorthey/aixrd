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