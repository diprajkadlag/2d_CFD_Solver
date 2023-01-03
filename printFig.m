function [] = printFig(phi,P,n,nOut,t)
if mod(n,nOut) == 0
    
    subplot(2,1,1)
    imagesc(phi);
    set(gca,'YDir','normal');
    colorbar;
    colormap(jet);
    title('phi');
     subplot(2,1,2)
    imagesc(P);
    set(gca,'YDir','normal');
    colorbar;
    colormap(jet);
    title('Pressure');
   
end

fprintf('Iteration No.: %d,  time: %d \n',n,t);
pause(0.0001)
end

   