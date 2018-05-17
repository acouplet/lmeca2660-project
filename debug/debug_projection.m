function debug_projection(Nx, dt)
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
    caxis([-1e-2 1e-2]);
    T = fread(fileID,[Nx+2,1.5*Nx+2],'double')';
    uspeedfile = fopen(sprintf('data/u_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    vspeedfile = fopen(sprintf('data/v_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    udebug = fopen(sprintf('data/ustar_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    vdebug = fopen(sprintf('data/vstar_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    u = fread(uspeedfile,[Nx+2,1.5*Nx+1],'double')';
    v = fread(vspeedfile,[Nx+1,1.5*Nx+2],'double')';
    us = fread(udebug,[Nx+2,1.5*Nx+1],'double')';
    vs = fread(vdebug,[Nx+1,1.5*Nx+2],'double')';
    fclose(uspeedfile);
    fclose(vspeedfile);
    fclose(udebug);
    fclose(vdebug);
    %quiver(X, Y, [u;zeros(1,Nx+2)],[v,zeros(Nx*1.5+2,1)], 'r');
    %hold on;
    %quiver(X, Y, [us;zeros(1,Nx+2)],[vs,zeros(Nx*1.5+2,1)], 'g');
    %hold off;
    %colorbar();
    subplot(2,2,1);
    contourf(X,Y,[u;zeros(1,Nx+2)],'LineStyle','none');
    subplot(2,2,2);
    contourf(X,Y,[us;zeros(1,Nx+2)],'LineStyle','none');
    subplot(2,2,3);
    contourf(X,Y,[v,zeros(Nx*1.5+2,1)],'LineStyle','none');
    subplot(2,2,4);
    contourf(X,Y,[vs,zeros(Nx*1.5+2,1)],'LineStyle','none');
    title(sprintf('k = %d',k));
    fclose(fileID);
    k = k + 50;
    sleep(1);
end


close all;

end
