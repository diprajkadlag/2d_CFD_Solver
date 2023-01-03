function [flxConXx,flxConYy] = calcFluxConUDS(phi,Us,Vs)

global Ifim  Ila  Jfim  Jla  Ifi Ilap Jfi Jlap dx



flxConXx = 0.5*dx^2*(+ max(0,Us(Ifim:Ila,:)).*phi(Ifim:Ila,:) + min(0,Us(Ifim:Ila,:)).*phi(Ifi:Ilap,:));
flxConYy = 0.5*dx^2*(+ max(0,Vs(:,Jfim:Jla)).*phi(:,Jfim:Jla) + min(0,Vs(:,Jfim:Jla)).*phi(:,Jfi:Jlap));

end