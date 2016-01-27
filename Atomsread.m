% mldread reads output in Molden format from Molpro
% Author: Thomas Northey % Date: 20/5/13
% Description: Generates inputs from 'mldfile', a Molden file (usually output from Molpro)

function [Atoms] = Atomsread(mldfile)

au2ang = 0.52917721092d0; % convert Ang coordinates to au at the end

Atoms=[];

disp(['Atomsread: Read from file ' mldfile])
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

tline = fgetl(fid);  
% Atoms:
while ischar(tline)    
    p = textscan(tline,'%s %f %f %f %f %f');   
    Atoms = [Atoms;cell2mat(p(2)) cell2mat(p(3)) ... 
        cell2mat(p(4))/au2ang...  % Atom Number, no of electrons, x, y, z (x,y,z in Bohr)
        cell2mat(p(5))/au2ang...
        cell2mat(p(6))/au2ang]; 
    br = strfind(tline,'[GTO]');           % Break on line containing '[GTO]'
    if ~isempty(br)
        break
    end 
    tline = fgetl(fid);                    % next line
end 
nelec=sum(Atoms(:,2));                     % Number of electrons
natom=max(Atoms(:,1));                     % Number of atoms

return  % END FUNCTION mldread
