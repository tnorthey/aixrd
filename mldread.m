% mldread reads output in Molden format from Molpro
% Author: Thomas Northey % Date: 20/5/13
% Description: Generates inputs from 'mldfile', a Molden file (usually output from Molpro)

function [bi,M,ga,ci,l,m,n,xx,yy,zz,ppmo,moocc,primlist,Atoms] = mldread(mldfile)

    au2ang = 0.52917721092d0; % convert Ang coordinates to au at the end
    
    cc=0; nmo=0; gcount=0;
    SS=[];PP=[];DD=[];FF=[];MM=[];Atoms=[];An=[];E=[];H=[];

    disp(['mldread: Read from file ' mldfile])
    fid = fopen(mldfile, 'r');                 % Open file for reading
    tline = fgetl(fid);                        % Read 1st line from file, removing newline characters
                                               % Go to line containing '[Atoms]'
    while ischar(tline)                       
        br = strfind(tline,'[Atoms]');          
    if ~isempty(br)                                 
        break                              
    end                                    
    tline = fgetl(fid);                    
    end

    % Atoms:
    while ischar(tline)    
        p = textscan(tline,'%*s %f %f %f %f %f');   
        Atoms = [Atoms;cell2mat(p(1)) cell2mat(p(2)) cell2mat(p(3))...
                 cell2mat(p(4)) cell2mat(p(5))];    % Atom Number, no of electrons, x, y, z        
        br = strfind(tline,'[GTO]');           % Break on line containing '[GTO]'
        if ~isempty(br)
            break
        end  
        tline = fgetl(fid);                    % next line
    end 
    nelec=sum(Atoms(:,2));                     % Number of electrons
    natom=max(Atoms(:,1));                     % Number of atoms

    % GTOs:
    while ischar(tline)                        % While loop
                                               % find atom number
        an = sscanf(tline,'  %d 0');
        if ~isempty(an)
            An=an;                                 % Atom number
        end
        % s-contractions
        ss = strfind(tline,' s ');             % Find line containing 
        if ~isempty(ss)                        % if ss is not empty (will be if does not contain ' s ')
            A = sscanf(tline,' s    %d');      % finds A, the number of primitives/contraction
            cc=cc+1;                           % counts from 1 to number of contractions
            S = zeros(1,7);                    % assign zero matrices
            for i = 1:A                        % loop A times
            gcount=gcount+1;               % counts primitives 
            tline = fgetl(fid);            % next line
            p = textscan(tline,'%f %f');   % read two floating points 
            S(i,:) = [cc An cell2mat(p(1)) cell2mat(p(2)) 0 0 0];   
            end
            SS=[SS;S];                         % write out to matrix 
        end         
        % p-contractions
        pp = strfind(tline,' p ');             % Find line containing 
        if ~isempty(pp)                        % if pp is not empty (will be if does not contain ' p ') 
            B = sscanf(tline,' p    %d');      % finds B, the number of primitives/contraction
            cc=cc+1;                           % counts from 1 to number of contractions
            Px=[];Py=[];Pz=[];                 % assign zero matrices        
            for i = 1:B                        % loop B times
                gcount=gcount+1;               % counts primitives 
                tline = fgetl(fid);            % next line
                q = textscan(tline,'%f %f');   % read two floating points 
                Px(i,:) =   [cc An cell2mat(q(1)) cell2mat(q(2)) 1 0 0];   
                Py(i,:) = [cc+1 An cell2mat(q(1)) cell2mat(q(2)) 0 1 0];
                Pz(i,:) = [cc+2 An cell2mat(q(1)) cell2mat(q(2)) 0 0 1];
            end
            PP=[PP;Px;Py;Pz];                  % write out to matrix 3 times (px, py, pz)
            cc=cc+2;                           % increase count appropriately
        end     
        % d-contractions
        dd = strfind(tline,' d ');             % Find line containing 
        if ~isempty(dd)                        % if dd is not empty (will be if does not contain ' d ')                
            C = sscanf(tline,' d    %d');      % finds C, the number of primitives/contraction
            cc=cc+1;                           % counts from 1 to number of contractions
            Dxx=[];Dyy=[];Dzz=[];              % assign zero matrices
            Dxy=[];Dxz=[];Dyz=[];
            for i = 1:C                        % loop C times
                gcount=gcount+1;               % counts primitives 
                tline = fgetl(fid);            % next line
                r = textscan(tline,'%f %f');   % read two floating points 
                Dxx(i,:) =   [cc An cell2mat(r(1)) cell2mat(r(2)) 2 0 0];      
                Dyy(i,:) = [cc+1 An cell2mat(r(1)) cell2mat(r(2)) 0 2 0];
                Dzz(i,:) = [cc+2 An cell2mat(r(1)) cell2mat(r(2)) 0 0 2];
                Dxy(i,:) = [cc+3 An cell2mat(r(1)) cell2mat(r(2)) 1 1 0];
                Dxz(i,:) = [cc+4 An cell2mat(r(1)) cell2mat(r(2)) 1 0 1];
                Dyz(i,:) = [cc+5 An cell2mat(r(1)) cell2mat(r(2)) 0 1 1];
            end
            DD=[DD;Dxx;Dyy;Dzz;Dxy;Dxz;Dyz];   % write out to matrix 6 times (d_ij)
            cc=cc+5;                           % increase count appropriately
        end     
        % f-contractions
        ff = strfind(tline,' f ');             % Find line containing 
        if ~isempty(ff)                        % if ff is not empty (will be if does not contain ' f ')        
            D = sscanf(tline,' f    %d');      % finds E, the number of primitives/contraction
            cc=cc+1;                           % counts from 1 to number of contractions
            Fxxx=[];Fyyy=[];Fzzz=[];           % assign zero matrices       
            Fxyy=[];Fxxy=[];Fxxz=[];
            Fxzz=[];Fyzz=[];Fyyz=[];Fxyz=[];
            for i = 1:D                        % loop E times
                gcount=gcount+1;               % counts primitives 
                tline = fgetl(fid);            % next line
                s = textscan(tline,'%f %f');   % read two floating points 
                Fxxx(i,:) =   [cc An cell2mat(s(1)) cell2mat(s(2)) 3 0 0];   
                Fyyy(i,:) = [cc+1 An cell2mat(s(1)) cell2mat(s(2)) 0 3 0];
                Fzzz(i,:) = [cc+2 An cell2mat(s(1)) cell2mat(s(2)) 0 0 3];
                Fxyy(i,:) = [cc+3 An cell2mat(s(1)) cell2mat(s(2)) 1 2 0];
                Fxxy(i,:) = [cc+4 An cell2mat(s(1)) cell2mat(s(2)) 2 1 0];
                Fxxz(i,:) = [cc+5 An cell2mat(s(1)) cell2mat(s(2)) 2 0 1];
                Fxzz(i,:) = [cc+6 An cell2mat(s(1)) cell2mat(s(2)) 1 0 2];
                Fyzz(i,:) = [cc+7 An cell2mat(s(1)) cell2mat(s(2)) 0 1 2];
                Fyyz(i,:) = [cc+8 An cell2mat(s(1)) cell2mat(s(2)) 0 2 1];
                Fxyz(i,:) = [cc+9 An cell2mat(s(1)) cell2mat(s(2)) 1 1 1];            
            end
            FF=[FF;Fxxx;Fyyy;Fzzz;Fxyy;Fxxy;Fxxz;Fxzz;Fyzz;Fyyz;Fxyz];  % write out to matrix 10 times (f_ijk)
            cc=cc+9;                           % increase count appropriately                    
        end               
        % Break on line containing '[MO]'    
        br = strfind(tline,'[MO]');           
        if ~isempty(br)
            break
        end    
        tline = fgetl(fid);        % Read next line    
    end     
    T=sortrows([SS;PP;DD;FF],1);   % T matrix
    ngc=max(T(:,1));               % Number of contractions
    ppmo=size(T,1);                % primitives per MO

    % MO coefficients & Occupancy numbers:
    while ischar(tline)   
        mos = strfind(tline,' Occup=');             % Find line containing 
        if ~isempty(mos)                            % if mos is not empty         
            bi = sscanf(tline,' Occup=    %f');     % finds bi (=0,1,2), no. of electrons per MO
            if bi>0
            nmo=nmo+1;                              % counts from 1 to number of MOs
            M = zeros(1,4);                         % assign zero matrix
            for i = 1:ngc                           % loop nmo times 
                tline = fgetl(fid);                 % next line
                s = textscan(tline,'%f %f');        % read MO number, MO coeff. 
                M(i,:) = [nmo cell2mat(s(1)) bi cell2mat(s(2))];  %
            end
            MM=[MM;M];            % write out to matrix
            moocc=max(nmo);
            end
        end 
        tline = fgetl(fid);       % Read next line
    end
    ngp=nmo*ppmo;                  % total no. of primitives
    for i=1:natom                  % combine 'T' and 'Atoms' matrices appropriately
        for j=1:ppmo
            if Atoms(i,1)==T(j,2)
                E=[E;T(j,:) Atoms(i,:)];    % makes matrix E with size ppmo
            end
        end
    end
    for k=1:nmo                    % makes H matrix like (some columns are the same):
        for i=(k-1)*ngc+1:k*ngc        % [MO# contraction# bi M contraction# atom# ga ci l m n...
            for j=1:ppmo               % atom# electrons/atom x y z], size nmo*ppmo = total no. of primitives.
                if MM(i,2)==E(j,1)
                    H=[H;MM(i,:) E(j,:)];
                end
            end
        end
    end
    % output list of atom indices associated with each primitive
    ntmp = length(H(1,:));
    primlist = H(:,ntmp-11); % atom index for each primitive, from
                             % T(j,2) via H-matrix length(H(:,1)) = ngp
    
    % Desired Variables
    bi=H(:,3);M=H(:,4);ga=H(:,7);ci=H(:,8);
    l=H(:,9);m=H(:,10);n=H(:,11);
    xx=H(:,14)/au2ang;yy=H(:,15)/au2ang;zz=H(:,16)/au2ang;
    
    % save('inputs')      % save to .mat file
    fclose(fid);
    disp('mldread: inputs generated')

end  % END FUNCTION mldread
