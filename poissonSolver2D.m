function [P,e,nIt] = poissonSolver2D(divP,P)

%  Solve (approximately) the Poisson equation for the pressure correction.
global dx Ifi Jfi Ila Jla Ifim Jfim Ilam Jlam Ifip Jfip Ilap Jlap dt

RHS         = divP/dt;     % RHS of poission equation
nIt         = 0;            % Initialising Number of iterations
e           = 100;           % Initialising error 
    


% max error and iterations
eMax      = 1e-3;
nItMax    = 30;

%% implementation of jacobi method

while(nIt<nItMax)&&(e>eMax)
    
    %% Increment iteration counter
    nIt         = nIt+1;
    %% pressure field of last iteration step
    P_old       = P(:,:);

    %% new pressure with Jacobi point iteration
    P(Ifi:Ila,Jfi:Jla)= 0.25...
               *(P_old(Ifim:Ilam,Jfi:Jla) + P_old(Ifip:Ilap,Jfi:Jla)...
                +P_old(Ifi:Ila,Jfip:Jlap) + P_old(Ifi:Ila,Jfim:Jlam)...
                   -RHS(Ifi:Ila,Jfi:Jla)*dx^2);

    %% boundary conditions 
    %Set zero gradient BC at the inlet and fixed value of zero at all
    P(:,Jfim)   = P(:,Jfi); % zero pressure gradient
  
    P(Ifim,:)   = 0;
    P(Ilap,:)   = 0;
    P(:,Jlap)   = 0;

    %% calculating max error
    e     = max(abs((P(Ifi:Ila,Jfi:Jla)-P_old(Ifi:Ila,Jfi:Jla))),[],'all');
        
       
end

fprintf('Pressure correction:eps= %d,  Pressure iterations= %d ', e, nIt);







    


