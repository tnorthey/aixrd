function [I,J] = find_significant_gtos3(M,c,ga,xx,yy,zz,ppmo,mo1,mo2,cutoff)
% T. Northey, 11/6/14
% Finds indices of GTO products with pre-coefficents > cutoff.
% This speeds up calculation as less GTOs to work with.
% 17/7/14: fixed bug from find_significant_gtos,
% re-titled as find_significant_gtos2
%=======================================
mos=mo2-mo1+1;
pre=zeros(mos*ppmo,mos*ppmo);
for mo=mo1:mo2                             % loop through MOs
    for i = (mo-1)*ppmo+1:mo*ppmo          % loop through GTOs
        for j = i:mo*ppmo                  % for j>=i
           
            % GTO overlap,
            gm=ga(i)+ga(j);
            K=exp(-ga(i)*ga(j)*gm^-1*((xx(i)-xx(j))^2+(yy(i)-yy(j))^2+(zz(i)-zz(j))^2));
            % pre-coeffs,
            pre(i,j)=K*M(i)*M(j)*c(i)*c(j);

        end
    end
end
%========================================
[I,J]=find(abs(pre)>cutoff);
% function ends
return

