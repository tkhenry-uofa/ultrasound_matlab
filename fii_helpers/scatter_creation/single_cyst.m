
function [positions, amp] = single_cyst (N,x_range, y_range, z_range, radius)


x_size = x_range(2)-x_range(1);   %  Width of phantom [mm]
y_size = y_range(2)-y_range(1);   %  Transverse width of phantom [mm]
z_size = z_range(2)-z_range(1);   %  Height of phantom [mm]


%  Create the general scatterers

x = rand (N,1)*x_size + x_range(1);
y = rand (N,1)*y_size + y_range(1);
z = rand (N,1)*z_size + z_range(1);


% Generate aplitudes with a Gaussian distribution

amp=randn(N,1);

%cyst
xc=x_size/2+x_range(1);  
yc=y_size/2+y_range(1);     %  Place of cyst [mm]
zc=z_size/2+z_range(1);  

inside = ( ((x-xc).^2 + (z-zc).^2) < radius^2);
% amp = amp .* (1-inside); 

x = x(not(inside));
y = y(not(inside));
z = z(not(inside));
amp = amp(not(inside));

%  Return the variables

positions=[x y z];



