function [C,Q,Qe,Qi,Fq,q] = setupq(Nq,wl,dim)
% setupq, setup q-vector grid 


Nq=2*ceil(Nq/2)-1;  % want Nq to be odd so q=0 is in linspace(-qmax,qmax,Nq)
                    % because q=0 is important (central maximum)

k0=2*pi/wl;         % incoming vector magnitude |k|
qmax=2*k0;          % q = k0 - k, where |k0|=|k| (elastic), so qmax=2|k|
if Nq==1
    q=0;
else
    q=linspace(0,qmax,Nq);  % define momentum transfer vector q
end

% coordinate choice (dim='x','y','z','xy','xz','yz','xyz','sph')
[qx,qy,qz,Fq,C] = choosecoord2(dim,q);

% if size(C,1)==3
%     r=C{1};
%     th=C{2};
%     ph=C{3};
% end

Q=cell(3,1);  % make convenient cell arrays
Qe=Q; Qi=Q;
Q{1}=qx; Q{2}=qy; Q{3}=qz;
for k=1:3
    Qe{k}=-.25*Q{k}.^2;
    Qi{k}=.5i*Q{k};
end

return

