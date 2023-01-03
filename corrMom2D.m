function [rhoU,rhoV] = corrMom2D(rhoU,rhoV,rhoUP,rhoVP,P)

global dx dt Ifi Jfi Ila Jla Ifip Jfip Ilap Jlap Ifim Jfim Jlam Ilam 

%% correction of rhoU

rhoU(Ifi:Ila,Jfi:Jla) = rhoUP(Ifi:Ila,Jfi:Jla)...
                        -dt*(P(Ifip:Ilap,Jfi:Jla)-P(Ifim:Ilam,Jfi:Jla))...
                        /(2.0 * dx);

%% correction of rhoV

rhoV(Ifi:Ila,Jfi:Jla) = rhoVP(Ifi:Ila,Jfi:Jla)...
                        -dt*(P(Ifi:Ila,Jfip:Jlap)-P(Ifi:Ila,Jfim:Jlam))...
                        /(2.0 * dx);


end

