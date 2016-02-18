function[Iq]=iamcalc2(Atoms,q)
% Approx. form factor equation from: 
% http://lamp.tu-graz.ac.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
% Uses Debye approx for orientational averaged diffraction pattern

% IAM rotationally-averaged diffraction intensity.
% au2ang = 0.52917721092d0;

Nq=length(q);        % length of q
Nat=size(Atoms,1);   % number of atoms
ZZ=Atoms(:,1);       % atomic number
X=Atoms(:,2); 
Y=Atoms(:,3);
Z=Atoms(:,4);

Iq=zeros(1,Nq);
for j1=1:Nat
    for j2=1:Nat           
        rij=sqrt((X(j1)-X(j2))^2+(Y(j1)-Y(j2))^2+(Z(j1)-Z(j2))^2);
        for i=1:Nq
            brac=q(i)*rij;
            if brac>1e-9
               Iq(i) = Iq(i) + formfact(ZZ(j1),q(i)).*formfact(ZZ(j2),q(i)).*sin(brac)./brac;
            else
               Iq(i) = Iq(i) + formfact(ZZ(j1),q(i)).*formfact(ZZ(j2),q(i));
            end    
        end
    end
end

return