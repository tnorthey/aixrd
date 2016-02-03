function[Fq,Atoms,C]=AIXRD2016_calcFq_iam(mldfile,Nq,wl,dim)

% Setup q-array,
[C,Q,~,~,Fq,q] = setupq(Nq,wl,dim);

% Read Atom information from molden file,
Atoms = Atomsread(mldfile);

Nat=Atoms(end,1);   % no. of atoms
R0=Atoms(:,3:5);    % nuclear geometry (Bohr)
an=Atoms(:,2);      % atomic numbers

%========================================
% calculate Fq

% Independent atom model method,
disp('Using IAM method...')
au2ang = 0.52917721092d0; % convert Ang coordinates to au
% total 3D form-factor, f(q),
if strfind(dim,'det')
    Fq = iam_detector(Atoms,dim,q);
else
    for i = 1:Nat           % loop through atoms         
        % Works with dim = {x, y, z, sph} only
        [fq,fq3]=formfact3(an(i),q/au2ang);
        phs=exp(1i*(Q{1}*R0(i,1)+Q{2}*R0(i,2)+Q{3}*R0(i,3)));
        if strcmp(dim,'x') || strcmp(dim,'y') || strcmp(dim,'z')
            Fq = Fq + fq'.*phs;    
        elseif strcmp(dim,'sph')
            Fq = Fq + fq3.*phs;   
        end
    end
end
    
%========================================
% For checking purposes,
maxFq=max(max(max(abs(Fq)))); % fmax = f(q=0) = Nelec
Nelec=maxFq;
Nelec0=sum(Atoms(:,2));       % Nelec from molden file
disp(strcat('True Nelec=',num2str(Nelec0)));
disp(strcat('Integrated Nelec=',num2str(Nelec))); % to check

return
