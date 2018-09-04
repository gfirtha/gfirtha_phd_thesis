clear
close all;

x0 = 1;
y0 = 0;
z0 = -1;

xs = 0;
ys = -2;
zs = 0;

x = 2; 
y = 2;
z = -2;

r1 = sqrt((x0-xs)^2+(y0-ys)^2+(z0-zs)^2);
r2 = sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2);
r3 = sqrt((x-xs)^2+(y-ys)^2+(z-zs)^2);
  
H1full = 1/r1*[    1-(x0-xs)^2/r1^2,     -(x0-xs)*(y0-ys)/r1^2,   -(x0-xs)*(z0-zs)/r1^2;
               -(x0-xs)*(y0-ys)/r1^2,     1-(y0-ys)^2/r1^2,  -(y0-ys)*(z0-zs)/r1^2;    
        -(x0-xs)*(z0-zs)/r1^2 ,       -(y0-ys)*(z0-zs)/r1^2,  1-(z0-zs)^2/r1^2 ];
kx1 = (x0-xs)/r1;
ky1 = (y0-ys)/r1;
kz1 = (z0-zs)/r1;

H2full = 1/r2*[    1-(x-x0)^2/r2^2,     -(x-x0)*(y-y0)/r2^2,   -(x-x0)*(z-z0)/r2^2;
               -(x-x0)*(y-y0)/r2^2,     1-(y-y0)^2/r2^2,  -(y-y0)*(z-z0)/r2^2;    
        -(x-x0)*(z-z0)/r2^2 ,       -(y-y0)*(z-z0)/r2^2,  1-(z-z0)^2/r2^2 ];
kx2 = (x-x0)/r2;
ky2 = (y-y0)/r2;
kz2 = (z-z0)/r2;

H3full = 1/r3*[    1-(x-xs)^2/r3^2,     -(x-xs)*(y-ys)/r3^2,   -(x-xs)*(z-zs)/r3^2;
               -(x-xs)*(y-ys)/r3^2,     1-(y-ys)^2/r3^2,  -(y-ys)*(z-zs)/r3^2;    
        -(x-xs)*(z-zs)/r3^2 ,       -(y-ys)*(z-zs)/r3^2,  1-(z-zs)^2/r3^2 ];

H = 1/r1*[    1-(x0-xs)^2/r1^2, -(x0-xs)*(z0-zs)/r1^2;    
        -(x0-xs)*(z0-zs)/r1^2 ,  1-(z0-zs)^2/r1^2 ];

[U1,V1] = eig(H1full);
[U2,V2] = eig(H2full);
[U3,V3] = eig(H3full);
[U,V] = eig(H1full+H2full);