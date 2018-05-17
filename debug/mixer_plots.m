Nx = 100;
Ny = Nx * 3 / 2;

L = 1;

h = L/Nx;

H = L * 3/2;

D = 3/5 * L;

a = D/2;

xc = L/2;
yc = L;

xhi = zeros(Nx, Ny);

for i = 1:Nx
  for j = 1:Ny
    x = h*i;
    y = h*j;
    theta = atan2( y - yc , x - xc );
    d = sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) );
    if ( d <= (a * cos(3 * theta))  )
      xhi(i,j) = 1;
    end
    %xhi(i,j) = ???
  end
end

% theta = linspace(0,2*pi,100);
% r = a * cos(3*theta);
%plot( L/2 + cos(theta).*r, L/2 + sin(theta).*r );
%hold on;
%plot( L/2 + cos(theta).*r/2, L/2 + sin(theta).*r/2 );

spy(xhi')