function plots_octave(Nx, dt)
close all;

Nx = 50;
X = linspace(0,1,Nx+2)
Y = linspace(1.5,0,1.5*Nx+2)
k = 0;
dt = 100;

while(1)
    fileID = fopen(sprintf('data/T_Nx%d_dt%d_iter%d.bin',Nx,dt,k))
    if(fileID == -1) break; end
    caxis manual;
    T = fread(fileID,[Nx+2,1.5*Nx+2],'double')';
    size(T)
    uspeedfile = fopen(sprintf('data/u_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    vspeedfile = fopen(sprintf('data/v_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    u = fread(uspeedfile,[Nx+2,1.5*Nx+1],'double')';
    v = fread(vspeedfile,[Nx+1,1.5*Nx+2],'double')';
    fclose(uspeedfile);
    fclose(vspeedfile);
    contourf(X,Y,T,'LineStyle','none')
    caxis([-1e-2 1e-2]);
    hold on;
    quiver(X, Y, [u;zeros(1,Nx+2)],[v,zeros(Nx*1.5+2,1)], 'k')
    hold off;
    colorbar();
    title(sprintf('k = %d',k));
    fclose(fileID);
    k = k + 50;
    sleep(1);
end


close all;

end
