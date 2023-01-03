function [Us,Vs]= mom2vel(rhoU,rhoV)
%mom2vel converts momentum into surface velocities

%calculates the velocity fields Us(:,:) and Vs(:,:) on the cell surfaces 
%from the momentum fields rhoU(:,:) and rhoV(:,:) that are stored
%on the cell centre and from a constant density rho

%Us(east)=(Uc(i+1,j)+Uc(i,j))/2

global rho Ifim Ifi Ila Ilap Jfim Jfi Jla Jlap 

%% interpolation of surface velocities
    Us  = (0.5/rho)*(rhoU(Ifim:Ila,:)+rhoU(Ifi:Ilap,:));
    Vs  = (0.5/rho)*(rhoV(:,Jfim:Jla)+rhoV(:,Jfi:Jlap));

end