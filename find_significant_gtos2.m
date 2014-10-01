function [I,J] = find_significant_gtos2(M,ci,ppmo,moocc,cutoff)
% T. Northey, 11/6/14
% Finds indices of GTO products with pre-coefficents > cutoff.
% This speeds up calculation as less GTOs to work with.
% 17/7/14: fixed bug from find_significant_gtos,
% re-titled as find_significant_gtos2
%=======================================
mo1=1; mo2=moocc;
pre=zeros(moocc*ppmo,moocc*ppmo);
for mo=mo1:mo2                             % loop through MOs
    for i = (mo-1)*ppmo+1:mo*ppmo          % loop through GTOs
        for j = (mo-1)*ppmo+1:mo*ppmo 
            pre(i,j)=M(i)*M(j)*ci(i)*ci(j);   
        end
    end
end
%========================================
[I,J]=find(abs(pre)>cutoff);
% function ends
return

