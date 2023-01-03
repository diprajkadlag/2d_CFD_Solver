function [rhoU,rhoV,phi,P] = initialisation()
%initialisation : initialises the fields rhoU, rhoV, phi ,P

global  rho Ilap Jlap velo wjet Imid 


%% Initialising fields with zeros
rhoU          = zeros(Ilap,Jlap); 
rhoV          = zeros(Ilap,Jlap);
phi           = zeros(Ilap,Jlap);
P             = zeros(Ilap,Jlap);
%% Defining momentum field 
 
 rhoV(int8(Imid-wjet/2-1):int8(Imid+wjet/2+1), :)  = rho*velo/2.0; 
 rhoV(int8(Imid-wjet/2):int8(Imid+wjet/2), :)      = rho*velo;
 
%% Defining scalar field 
 phi(int8(Imid-wjet/2-1):int8(Imid+wjet/2+1), :)   = 0.5*rho;
 phi(int8(Imid-wjet/2):int8(Imid+wjet/2), :)       = 1.0*rho;
 
end


