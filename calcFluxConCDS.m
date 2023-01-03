function [flxConX,flxConY] = calcFluxConCDS(phi,Us,Vs)
% calcFlxConCDS calculates the convective fluxes fluxConX and fluxConY of a conserved 
% quantity phi in x and y direction over the e and n surfaces of each cell.
% Using central differencing scheme (CDS)

global rho Ifim  Ila  Jfim  Jla  Ifi Ilap Jfi Jlap
%% convective fluxes in X and Y direction 
  flxConX   = 0.5*rho * (phi(Ifi:Ilap,:)+phi(Ifim:Ila,:)).* Us(Ifim:Ila,:);
  flxConY   = 0.5*rho * (phi(:,Jfi:Jlap)+phi(:,Jfim:Jla)).* Vs(:,Jfim:Jla);
 
end