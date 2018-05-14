function plots(Nx, dt, saveIter, mixer)
close all;

X = linspace(0,1,Nx+2)
Y = linspace(1.5,0,1.5*Nx+2)
k = 0;
figure('Position',[10 10 600 600]);
v = VideoWriter('test.avi');
v.FrameRate = 10;
open(v)
while(1)
    fileID = fopen(sprintf('data/T_Nx%d_dt%d_iter%d_mixing%d.bin',Nx,dt,k,mixer))
    if(fileID == -1) break; end
    T = fread(fileID,[Nx+2,1.5*Nx+2],'double')';
    contourf(X,Y,T,'LineStyle','none');
    colorbar();
    title(sprintf('k = %d',k));
    writeVideo(v,getframe(1)); 
    fclose(fileID);
    k = k + saveIter;
end


close(v)
close all;

end
