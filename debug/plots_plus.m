function plots_plus(Nx, dt)
close all;

Ny = 1.5*Nx;
h = (1/1.5)/Nx;
XT = linspace(-(h/2),1/1.5 + h/2,Nx+2);
YT = linspace(1+h/2,-h/2,1.5*Nx+2);
k = 0;
%figure('Position',[10 10 600 600]);
v = VideoWriter('test.avi');
v.FrameRate = 10;
open(v)
while(1)
    fileID = fopen(sprintf('data/T_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    if(fileID == -1) break; end
    T = fread(fileID,[Nx+2,1.5*Nx+2],'double')';
    %fclose(fileID);
    %fileID = fopen(sprintf('data/u_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    %if(fileID == -1) break; end
    %ux = fread(fileID,[Nx+2,1.5*Nx+1],'double')';
    %fclose(fileID);
    %fileID = fopen(sprintf('data/mix_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    %if(fileID == -1) break; end
    %mixer = fread(fileID,[Nx,1.5*Nx],'double')';
    %fclose(fileID);
    %fileID = fopen(sprintf('data/v_Nx%d_dt%d_iter%d.bin',Nx,dt,k));
    %if(fileID == -1) break; end
    %caxis manual;
    %uy = fread(fileID,[Nx+1,1.5*Nx+2],'double')';
    %fclose(fileID);
    
    pcolor(XT,YT,T);
    xlim([1.9/5/1.5,3.1/5/1.5]);
    ylim([1.9/5,3.1/5]);
    axis equal;
    %contourf(XT,YT,T,'LineStyle','none');
    %contourf(XT,YT,[ux;zeros(1,Nx+2)]);
    hold on;
    %contourf(X(2:Nx+1),Y(2:1.5*Nx+1),mixer,'LineStyle', 'none');% TOFIX: values between 0 and 1, way out of the range :(
    %quiver(X, Y, [ux;zeros(1,Nx+2)],[uy,zeros(Nx*1.5+2,1)], 'k');
    hold off;
    colorbar();
    caxis([0 10]);
    title(sprintf('k = %d',k));
    writeVideo(v,getframe(1)); 

    %pause();
    k = k + 1;
end


close(v)
close all;

end
