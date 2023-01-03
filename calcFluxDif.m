function [flxDifX,flxDifY] = calcFluxDif(phi)
% diffusive fluxes 'fluxDifX' and 'fluxDifY' of a conserved scalar phi in
% x and y direction over the e and n surfaces of each cell, depending on phi,
% D and dx.

global rho dx D Ifim Ifi Ila Ilap Jfim Jfi Jla Jlap

%% diffusive fluxes in x and y direction
 flxDifX    = (D*dx/rho)*(phi(Ifi:Ilap,:)-phi(Ifim:Ila,:));
 flxDifY    = (D*dx/rho)*(phi(:,Jfi:Jlap)-phi(:,Jfim:Jla)); 

end
