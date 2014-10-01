function [ft] = ftswch2(ga1,ga2,L1,L2,x1,x2,s)
% Author: Thomas Northey % Date: 22/03/13

% DESCRIPTION:
% Function calculates FT products for all s..f-type GTOs
% given appropriate L=(l1,m1,n1),(l2,m2,n2) values, exponents (ga1,ga2), 
% s-vector (sx,sy,sz), and GTO position x=(x1,y1,z1),(x2,y2,z2)

% EDIT: 17/06/14. Slightly faster than ftswch. No pre-coefficients.

switch L1
    case 0
        switch L2
            case 0
                ft = ft00(ga1,ga2,x1,x2,s);
            case 1
                ft = ft10(ga2,ga1,x2,x1,s);
            case 2
                ft = ft20(ga2,ga1,x2,x1,s);
            case 3
                ft = ft30(ga2,ga1,x2,x1,s);
        end
    case 1
        switch L2
            case 0
                ft = ft10(ga1,ga2,x1,x2,s);
            case 1
                ft = ft11(ga1,ga2,x1,x2,s);
            case 2
                ft = ft21(ga2,ga1,x2,x1,s);
            case 3
                ft = ft31(ga2,ga1,x2,x1,s);
        end
    case 2
        switch L2
            case 0
                ft = ft20(ga1,ga2,x1,x2,s);
            case 1
                ft = ft21(ga1,ga2,x1,x2,s);
            case 2
                ft = ft22(ga1,ga2,x1,x2,s);
            case 3
                ft = ft32(ga2,ga1,x2,x1,s);  
        end   
    case 3
        switch L2
            case 0
                ft = ft30(ga1,ga2,x1,x2,s);
            case 1
                ft = ft31(ga1,ga2,x1,x2,s);
            case 2
                ft = ft32(ga1,ga2,x1,x2,s); 
            case 3
                ft = ft33(ga1,ga2,x1,x2,s);
        end
end
return

% FUNCTIONS:
% ss GTO product, calculates FT[g00], where g0=exp(-ga(x-x0)^2): 
function[FT00] = ft00(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
norm1 = (2*ga1*pi^-1)^.25;
norm2 = (2*ga2*pi^-1)^.25;
FT00 = K*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return

% ps GTO product, calculates FT[g10]:
function[FT10] = ft10(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
pre2 = (0.5i*s*ga3^-1+a-x1);
norm1 = 2*ga1^.5*(2*ga1*pi^-1)^.25;
norm2 = (2*ga2*pi^-1)^.25;
FT10 = K*pre2*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return

% pp GTO product, calculates FT[g11]:
function[FT11] = ft11(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
A = 0.5i*s*ga3^-1+a;
pre2 = (A^2+0.5*ga3^-1-(x1+x2)*A+x1*x2);
norm1 = 2*ga1^.5*(2*ga1*pi^-1)^.25;
norm2 = 2*ga2^.5*(2*ga2*pi^-1)^.25;
FT11 = K*pre2*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return

% ds GTO product, calculates FT[g20]:
function[FT20] = ft20(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
A = 0.5i*s*ga3^-1+a;
pre2 = (A^2+0.5*ga3^-1-2*x1*A+x1^2);
norm1 = (4*ga1*3^-.5)*(2*ga1*pi^-1)^.25;
norm2 = (2*ga2*pi^-1)^.25;
FT20 = K*pre2*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return

% dp GTO product, calculates FT[g21]:
function[FT21] = ft21(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
A = 0.5i*s*ga3^-1+a;
pre2 = (A^3-(2*x1+x2)*A^2+(x1^2+2*x1*x2+1.5*ga3^-1)*A-(2*x1+x2)*0.5*ga3^-1-x1^2*x2);
norm1 = (4*ga1*3^-.5)*(2*ga1*pi^-1)^.25;
norm2 = (32*ga2^3*pi^-1)^.25;
FT21 = K*pre2*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return

% dd GTO product, calculates FT[g22]:
function[FT22]=ft22(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
A = 0.5i*s*ga3^-1+a;
pre2 = (A^4-2*(x1+x2)*A^3+(3*ga3^-1+x1^2+x2^2+4*x1*x2)*A^2-(2*x1^2*x2+...
    2*x1*x2^2+3*(x1+x2)*ga3^-1)*A+0.5*(x1^2+x2^2+4*x1*x2)*ga3^-1+0.75*ga3^-2+x1^2*x2^2);
norm1 = (4*ga1*3^-.5)*(2*ga1*pi^-1)^.25;
norm2 = (4*ga2*3^-.5)*(2*ga2*pi^-1)^.25;
FT22 = K*pre2*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return

% fs GTO product, calculates FT[g30]:
function[FT30]=ft30(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
A = 0.5i*s*ga3^-1+a;
pre2 = (A^3-3*x1*A^2+(3*x1^2+1.5*ga3^-1)*A-1.5*x1*ga3^-1-x1^3);
norm1 = (8*ga1^1.5*5^-.5)*(2*ga1*pi^-1)^.25;
norm2 = (2*ga2*pi^-1)^.25;
FT30 = K*pre2*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return
 
% fp GTO product, calculates FT[g31]:
function[FT31]=ft31(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
A = 0.5i*s*ga3^-1+a;
pre2 = (A^4-(3*x1+x2)*A^3+(3*x1^2+3*x1*x2+3*ga3^-1)*A^2-(3*x1^2*x2+x1^3+...
    1.5*x2*ga3^-1+4.5*x1*ga3^-1)*A+(3*x1*x2+3*x1^2)*(2*ga3)^-1+0.75*ga3^-2+x1^3*x2);
norm1 = (8*ga1^1.5*5^-.5)*(2*ga1*pi^-1)^.25;
norm2 = (32*ga2^3*pi^-1)^.25;
FT31 = K*pre2*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return

% fd GTO product, calculates FT[g32]:
function[FT32]=ft32(ga1,ga2,x1,x2,s)
ga3 = ga1+ga2;
K = exp(-ga1*ga2*ga3^-1*(x1-x2)^2);
a = (ga1*x1+ga2*x2)*ga3^-1;
% A = 0.5i*s*ga3^-1+a;
pre2 = (32*ga3^5)^-1*...
    (32*a^5*ga3^5 + 1i*s^5 - 32*ga3^5*x1^3*x2^2 -...
    16*a^4*ga3^4*(-5i*s + 6*ga3*x1 + 4*ga3*x2) + 2*ga3*1i*s^3*(-10 + 1i*s*(3*x1 + 2*x2)) +...
    16*ga3^4*x1*(-3*x2^2 + 3*x1*x2*(-2 + 1i*s*x2) + x1^2*(-1 + 2i*s*x2)) + 16*a^3*ga3^3*...
    (-5*s^2 - 2*ga3*(-5 + 6i*s*x1 +4i*s*x2) + 2*ga3^2*(3*x1^2 + 6*x1*x2 + x2^2)) +...
    4i*ga3^2*s*(15 - 6i*s*(3*x1 + 2*x2) - s^2*(3*x1^2 + 6*x1*x2 +x2^2)) -...
    8*ga3^3*(-s^2*x1^3 - 3*x2*(-2 + 1i*s*x2) + 3i*s*x1^2*(-3 + 2i*s*x2) +...
    3*x1*(3 - 6i*s*x2 - s^2*x2^2)) - 8*a^2*ga3^2*(5i*s^3 +...
    4*ga3^3*x1*(x1^2 + 6*x1*x2 + 3*x2^2) + 6i*ga3*s*(-5 + 1i*s*(3*x1 + 2*x2)) -...
    6*ga3^2*(3i*s*x1^2 + x2*(-4 + 1i*s*x2) + 6*x1*(-1 + 1i*s*x2))) +...
    2*a*ga3*(5*s^4 + 16*ga3^4*x1^2*x2*(2*x1 + 3*x2) +...
    4*ga3*s^2*(-15 + 6i*s*x1 + 4i*s*x2) -...
    8*ga3^3*(2i*s*x1^3 - 3*x2^2 + 6*x1*x2*(-3 +1i*s*x2) + 3*x1^2*(-3 + 4i*s*x2))+...
    12*ga3^2*(5 - 4i*s*(3*x1 + 2*x2) - s^2*(3*x1^2 + 6*x1*x2 + x2^2))));
norm1 = (8*ga1^1.5*5^-.5)*(2*ga1*pi^-1)^.25;
norm2 = (4*ga2*3^-.5)*(2*ga2*pi^-1)^.25;
FT32 = K*pre2*norm1*norm2*(pi/ga3)^.5*exp(-0.25*s^2*ga3^-1)*(cos(s*a)+1i*sin(s*a));
return

% ff GTO product, calculates FT[g33]:
function[FT33]=ft33(ga1,ga2,x1,x2,s)
FT33=0;   % Algebra became too difficult.. f-f overlap ~ 0 in most cases.
return

% % gs GTO product, calculates FT[g40]:
% function[FT40]=ft40(ga1,ga2,x1,x2,s) 
% ga3 = ga1+ga2;
% a = (ga1*x1+ga2*x2)*ga3^-1;
% A = 0.5i*s*ga3^-1+a;
% pre = (A^4-4*x1*A^3+(3*ga3^-1+6*x1^2)*A^2-(4*x1^3+6*x1*ga3^-1)*A+3*x1^2*ga3^-1+x1^4+0.75*ga3^-2);
% FT40 = pre*ft00(ga1,ga2,x1,x2,s);
% return