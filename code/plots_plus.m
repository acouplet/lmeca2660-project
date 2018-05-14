function plots_plus(Nx, dt)
close all;

X = linspace(0,1,Nx+2);
Y = linspace(1.5,0,1.5*Nx+2);
k = 0;
figure('Position',[10 10 600 600]);
v = VideoWriter('test.avi');
v.FrameRate = 10;
open(v)
while(1)
    fileID = fopen(sprintf('data/T_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    if(fileID == -1) break; end
    T = fread(fileID,[Nx+2,1.5*Nx+2],'double')';
    fclose(fileID);
    fileID = fopen(sprintf('data/u_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    if(fileID == -1) break; end
    ux = fread(fileID,[Nx+2,1.5*Nx+1],'double')';
    fclose(fileID);
    fileID = fopen(sprintf('data/mix_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    if(fileID == -1) break; end
    mixer = fread(fileID,[Nx,1.5*Nx],'double')';
    fclose(fileID);
    fileID = fopen(sprintf('data/v_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    if(fileID == -1) break; end
    caxis manual;
    uy = fread(fileID,[Nx+1,1.5*Nx+2],'double')';
    fclose(fileID);
    
    contourf(X,Y,T,'LineStyle','none');
    hold on;
    contourf(X(2:Nx+1),Y(2:1.5*Nx+1),mixer,'LineStyle', 'none');
    quiver(X, Y, [ux;zeros(1,Nx+2)],[uy,zeros(Nx*1.5+2,1)], 'k');
    hold off;
    colorbar();
    caxis([-1e-2 1e-2]);
    title(sprintf('k = %d',k));
    writeVideo(v,getframe(1)); 

    k = k + 50;
end


close(v)
close all;

end
