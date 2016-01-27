function [qx,qy,qz,Fq,C] = choosecoord2(dim,q)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% defaults
Nq=length(q);
C=q;

r=0; th=0; ph=0;    

if strcmp(dim,'x')
    % 1d cartesian coordinate system    
    qx=q;qy=0;qz=0;
    Fq=zeros(1,Nq);
elseif strcmp(dim,'y')
    % 1d cartesian coordinate system    
    qx=0;qy=q;qz=0;
    Fq=zeros(1,Nq);
elseif strcmp(dim,'z')
    % 1d cartesian coordinate system    
    qx=0;qy=0;qz=q;
    Fq=zeros(1,Nq);
elseif strcmp(dim,'xy')
    % 2d cartesian coordinate system
    [qx,qy]=meshgrid(q);
    qz=0;
    Fq=zeros(Nq,Nq);    
elseif strcmp(dim,'xz')
    % 2d cartesian coordinate system
    [qx,qz]=meshgrid(q);
    qy=0;
    Fq=zeros(Nq,Nq); 
elseif strcmp(dim,'yz')
    % 2d cartesian coordinate system
    [qy,qz]=meshgrid(q);
    qx=0;
    Fq=zeros(Nq,Nq);     
elseif strcmp(dim,'xyz')
    % 3d cartesian coordinate system
    [qx,qy,qz]=meshgrid(q);
    Fq=zeros(Nq,Nq,Nq);
elseif strcmp(dim,'sph');
    Fq=zeros(Nq,Nq,Nq);
    %=======================================================================
    % spherical coordinate system
    nr=Nq;
    nth=Nq;
    nph=Nq;

    qmax=max(q);
    r=linspace(0,qmax,nr);       % radius
    th=linspace(0,pi,nth);       % theta
    ph=linspace(0,2*pi,nph);     % phi

    qx=zeros(nr,nth,nph);qy=qx;qz=qx;

    for k=2:nr                     % loop through spheres of non-zero radius r(kk)    
        for i=1:nph                % phi loop, note: skips 2*pi as f(0)=f(2*pi)
            for j=1:nth            % theta loop           
                % Create x,y,z as a function of spherical coords...           
                qx(k,j,i)=r(k)*sin(th(j))*cos(ph(i));
                qy(k,j,i)=r(k)*sin(th(j))*sin(ph(i));
                qz(k,j,i)=r(k)*cos(th(j));        
            end
        end  
    end
    C=cell(3,1); C{1}=r; C{2}=th; C{3}=ph; % coordinate matrix
    %=======================================================================
elseif strcmp(dim,'detx');
    % x-ray along x-axis
    Fq=zeros(Nq,Nq);
    %=======================================================================
    % spherical detector coordinate system
    nth=Nq;
    nph=Nq;

    qmax=max(q);
    r=qmax;         % radius
    th=linspace(0,pi,nth);         % theta
    ph=linspace(0,2*pi,nph);     % phi

    qx=zeros(nth,nph);qy=qx;qz=qx;
    
    for i=1:nph                % phi loop, note: skips 2*pi as f(0)=f(2*pi)
        for j=1:nth            % theta loop           
            % Create x,y,z as a function of spherical coords...           
            qy(j,i)=r*sin(th(j))*cos(ph(i));
            qz(j,i)=r*sin(th(j))*sin(ph(i));
            qx(j,i)=r*cos(th(j))-r;        
        end
    end  
    C=cell(3,1); C{1}=r; C{2}=th; C{3}=ph; % coordinate matrix
    %======================================================================= 
elseif strcmp(dim,'dety');
    % x-ray along y-axis
    Fq=zeros(Nq,Nq);
    %=======================================================================
    % spherical detector coordinate system
    nth=Nq;
    nph=Nq;

    qmax=max(q);
    r=qmax;         % radius
    th=linspace(0,pi,nth);         % theta
    ph=linspace(0,2*pi,nph);     % phi

    qx=zeros(nth,nph);qy=qx;qz=qx;
    
    for i=1:nph                % phi loop, note: skips 2*pi as f(0)=f(2*pi)
        for j=1:nth            % theta loop           
            % Create x,y,z as a function of spherical coords...           
            qz(j,i)=r*sin(th(j))*cos(ph(i));
            qx(j,i)=r*sin(th(j))*sin(ph(i));
            qy(j,i)=r*cos(th(j))-r;        
        end
    end 
    C=cell(3,1); C{1}=r; C{2}=th; C{3}=ph; % coordinate matrix
    %=======================================================================
elseif strcmp(dim,'detz');
    % x-ray along z-axis
    Fq=zeros(Nq,Nq);
    %=======================================================================
    % spherical detector coordinate system
    nth=Nq;
    nph=Nq;

    qmax=max(q);
    r=qmax;         % radius
    th=linspace(0,pi,nth);         % theta
    ph=linspace(0,2*pi,nph);     % phi

    qx=zeros(nth,nph);qy=qx;qz=qx;
    
    for i=1:nph                % phi loop, note: skips 2*pi as f(0)=f(2*pi)
        for j=1:nth            % theta loop           
            % Create x,y,z as a function of spherical coords...           
            qx(j,i)=r*sin(th(j))*cos(ph(i));
            qy(j,i)=r*sin(th(j))*sin(ph(i));
            qz(j,i)=r*cos(th(j))-r;        
        end
    end 
    C=cell(3,1); C{1}=r; C{2}=th; C{3}=ph; % coordinate matrix
    %======================================================================= 
elseif strcmp(dim,'detz0');
    % x-ray along z-axis
    Fq=zeros(Nq,Nq);
    %=======================================================================
    % spherical detector coordinate system
    nth=Nq;
    nph=Nq;

    qmax=max(q);
    r=qmax;         % radius
    th=linspace(0,pi,nth);         % theta
    ph=linspace(0,2*pi,nph);     % phi

    qx=zeros(nth,nph);qy=qx;qz=qx;
    
    for i=1:nph                % phi loop, note: skips 2*pi as f(0)=f(2*pi)
        for j=1:nth            % theta loop           
            % Create x,y,z as a function of spherical coords...           
            qx(j,i)=r*sin(th(j))*cos(ph(i));
            qy(j,i)=r*sin(th(j))*sin(ph(i));
            qz(j,i)=r*cos(th(j))-r;        
        end
    end 
    C=cell(3,1); C{1}=r; C{2}=th; C{3}=ph; % coordinate matrix
    %=======================================================================    
end

return

