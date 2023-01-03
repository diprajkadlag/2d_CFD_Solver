function [divg] = divMom2D(rhoUP,rhoVP)
global dx Ifi Ila Jfi Jla Ifip Ilap Ifim Ilam Jfip Jfim Jlam Jlap

%% divergence of the flow field, div = du/dx + dv/dy (2D)
divg(Ifi:Ila,Jfi:Jla) =(0.5/dx)...
                     *(rhoUP(Ifip:Ilap,Jfi:Jla)-rhoUP(Ifim:Ilam,Jfi:Jla)...
                       +rhoVP(Ifi:Ila,Jfip:Jlap)-rhoVP(Ifi:Ila,Jfim:Jlam) );

end