clc; clear;

%% Defining global variables
global Ima Jma dx rho D CFL nG Ifim Ifi Ifip Ilam Ila Ilap Jfim Jfi Jfip...
       Jlam Jla Jlap velo wjet dt Imid mu
  
%% parameters
Ima               = 128;            % number of cells vertically
Jma               = 256;            % number of cells horizontally
dx                = 5e-3;           % cell thickness dx=dy 
rho               = 1.2;            % density
mu                = 2e-5;           % dynamic viscocity
Sc                = 1.2;            % Schmidt-Number       
D                 = mu/(rho*Sc);    % diffusivity
CFL               = 0.7;           % Courant-Friedrich-Levy criterion CFL<1


%% Iteration and output parameters
MaxIt             = 10000;          % max no of iterations
n                 = 0; 
nOut              = 20;             % for printing o/p every 10 iterations
t                 = 0;

%% Defining convention
nG                 = 1;             % no of ghost cells
% vertical direction
Ifi                = nG+1;          % first cell in the numerical domain 
                                    % in vertical direction
Ifim               = nG;            % m stands for 'minus 1'
Ifip               = nG+2;          % p stands for 'plus 1'
Ila                = nG+Ima;        % last cell - vertically
Ilam               =(nG+Ima)-1;
Ilap               =(nG+Ima)+1;
Imid               = ceil(Ilap/2);  % mid cell - vertically 
% horizontal direction 
Jfi                = nG+1;
Jfim               = nG;
Jfip               = nG+2;
Jla                = nG+Jma;
Jlam               = (nG+Jma)-1;
Jlap               = (nG+Jma)+1;



%% initialisation of momemtum and scalar field
velo               = 1;            % velocity of the jet
wjet               = 24;           % jet width (2D)
[rhoU,rhoV,phi,P]  = initialisation();
rhoUP=rhoU;
rhoVP=rhoV;

%% velocity at cell surface from cell center 
[Us,Vs]           = mom2vel(rhoU,rhoV);

%% time-step width 
U                  = rhoU/rho;
V                  = rhoV/rho;
for i = Ifim:Ilap
    for j = Jfim:Jlap
        VeloAbs(i,j)= sqrt(U(i,j)^2 + V(i,j)^2);  %absolute velocity 
    end
end
VeloMx            = max(VeloAbs,[],'all');        % max absolute velocity
dt                = CFL*dx/VeloMx;                % dt=CFL*dx/max(|u,v,w|)

%% Iterations starts from here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i=1:MaxIt
    
%% convective fluxes of phi
[flxConX,flxConY] = calcFluxConUDS(phi,Us,Vs);
%% diffusive fluxes of phi
[flxDifX,flxDifY] = calcFluxDif(phi);
%% convective fluxes-rhoU
[flxConX_rhoU,flxConY_rhoU] = calcFluxConUDS(rhoU,Us,Vs);
%% diffusive fluxes-rhoU
[flxDifX_rhoU,flxDifY_rhoU] = calcFluxDif(rhoU*mu/D);
%% convective fluxes-rhoV
[flxConX_rhoV,flxConY_rhoV] = calcFluxConUDS(rhoV,Us,Vs);
%% diffusive fluxes-rhoV
[flxDifX_rhoV,flxDifY_rhoV] = calcFluxDif(rhoV*mu/D);


%% Transport of phi ***************************************************** 

%% phi in next timestep 
phi(Ifi:Ila,Jfi:Jla) = phi(Ifi:Ila,Jfi:Jla)...
                       +(dt/dx^3)*( (flxConX(Ifim:Ilam,Jfi:Jla)...
                                     - flxConX(Ifi:Ila,Jfi:Jla) )...
                                   +(flxConY(Ifi:Ila,Jfim:Jlam)...
                                     - flxConY(Ifi:Ila,Jfi:Jla) ) )...
                                     ...
                      -(dt/dx^3)*( (flxDifX(Ifim:Ilam,Jfi:Jla)...
                                    - flxDifX(Ifi:Ila,Jfi:Jla) )...
                                  +(flxDifY(Ifi:Ila,Jfim:Jlam)...
                                    - flxDifY(Ifi:Ila,Jfi:Jla) ) );                   
                        
%% Transport of rhoU ****************************************************

%% rhoU in the next step
rhoUP(Ifi:Ila,Jfi:Jla) = rhoU(Ifi:Ila,Jfi:Jla)...
                         +(dt/dx^3)*( (flxConX_rhoU(Ifim:Ilam,Jfi:Jla)...
                                       - flxConX_rhoU(Ifi:Ila,Jfi:Jla))...
                                     +(flxConY_rhoU(Ifi:Ila,Jfim:Jlam)...
                                       - flxConY_rhoU(Ifi:Ila,Jfi:Jla)))...
                                     ...
                        -(dt/dx^3)*( (flxDifX_rhoU(Ifim:Ilam,Jfi:Jla)...
                                      - flxDifX_rhoU(Ifi:Ila,Jfi:Jla))...
                                    +(flxDifY_rhoU(Ifi:Ila,Jfim:Jlam)...
                                      - flxDifY_rhoU(Ifi:Ila,Jfi:Jla)));
                          
%% Transport of rhoV ****************************************************

%% rhoV in next step
rhoVP(Ifi:Ila,Jfi:Jla) = rhoV(Ifi:Ila,Jfi:Jla)...
                        +(dt/dx^3)*( (flxConX_rhoV(Ifim:Ilam,Jfi:Jla)...
                                      -flxConX_rhoV(Ifi:Ila,Jfi:Jla))...
                                    +(flxConY_rhoV(Ifi:Ila,Jfim:Jlam)...
                                      - flxConY_rhoV(Ifi:Ila,Jfi:Jla)))...
                                     ...
                       -(dt/dx^3)*( (flxDifX_rhoV(Ifim:Ilam,Jfi:Jla)...
                                     - flxDifX_rhoV(Ifi:Ila,Jfi:Jla))...
                                   +(flxDifY_rhoV(Ifi:Ila,Jfim:Jlam)...
                                     - flxDifY_rhoV(Ifi:Ila,Jfi:Jla)));

%% projection method *****************************************************
% Calculate the pressure term using a projection method 
% Such that there is no divergence in the flow field!
%% divergence of flow-field
[divP]           = divMom2D(rhoUP,rhoVP);

%% pressure correction
% The pressure field can be computed fr[m the div in the velocity field!
[P,e,nIt]       = poissonSolver2D(divP,P);

%% momentum correction for the next step
[rhoU,rhoV]   = corrMom2D(rhoU,rhoV,rhoUP,rhoVP,P);
n                 = n + 1 ;
t                 = t + dt;

%% Boundary conditions

% phi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi(Ifim,:)   = phi(Ifi,:); %  TOP
phi(Ilap,:)   = phi(Ila,:); %  BOTTOM
phi(:,Jfim)   = 0; %  LEFT
phi(:,Jlap)   = phi(:,Jla); % von Neumann boundary condition RIGHT

phi(int8(Imid-wjet/2-1):int8(Imid+wjet/2+1), Jfim)   = 0.5*rho;
phi(int8(Imid-wjet/2):int8(Imid+wjet/2), Jfim)       = 1.0*rho;
 
%rhoU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoU(Ifim,:)          = 0;                       % tOP
rhoU(Ilap,:)          = 0;                       % BOTTOM
rhoU(:,Jfim)          = 0;
rhoU(:,Jfim)          = 0.5*(rand(Ilap,1)-0.5); %random noise at the inlet
                                                 %(vertical momentum)LEFT
rhoU(:,Jlap)          = rhoU(:,Jla);           % RIGHT

%rhov%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoV(Ifim,:)          = 0;                      % TOP
rhoV(Ilap,:)          = 0;                      % BOTTOM
rhoV(:,Jfim)          = 0;                      % LEFT
rhoV(:,Jlap)          = max(rhoV(:,Jla), 0);    % von Neumann BC

rhoV(int8(Imid-wjet/2-1):int8(Imid+wjet/2+1), Jfim)  = rho*velo/2; 
rhoV(int8(Imid-wjet/2):int8(Imid+wjet/2), Jfim)      = rho*velo;


%% velocity at cell surface from cell center 
[Us,Vs]           = mom2vel(rhoU,rhoV);

%% printing output after every 10 timesteps
printFig(phi,P,n,nOut,t);

end
set(gca, 'Ydir', 'normal');
