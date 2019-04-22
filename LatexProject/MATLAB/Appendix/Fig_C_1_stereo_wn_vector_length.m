clear
close all
addpath(genpath('../Files/arrow3'));
addpath(genpath('../Files'));
dx = 2.5e-3;
x = (-2:dx:2)';
y = (-1:dx:2.5)';

w = 0.8e3*2*pi;
c = 343.1;
k = w / c;
lambda = 2*pi/k;

R0 = 2.5;
fi0 = 30;

x1 = R0*[-sind(fi0) cosd(fi0)];
x2 = R0*[ sind(fi0) cosd(fi0)];

fi = 10;
C = ( tan(fi*pi/180)/tand(fi0) );
A1 = 1;
A2 = A1*(1-C)/(C+1);
%
Avec = [A1 A2];
A1 = A1/norm(Avec);  A2 = A2/norm(Avec);

[X,Y] = meshgrid(x,y);
r1 = sqrt((X-x1(1)).^2 + (Y-x1(2)).^2);
r2 = sqrt((X-x2(1)).^2 + (Y-x2(2)).^2);

r1px = (X-x1(1))./r1;                       r2px = (X-x2(1))./r2;
kx = (k*(r1px.*A1^2.*r2.^2+ r2px.*A2^2.*r1.^2 + A1*A2*r1.*r2.*(r1px+r2px).*cos(k*(r1-r2)))...-
    - A1*A2*(r1px.*r2-r2px.*r1).*sin(k*(r1-r2)))./((A1*r2).^2+(A2*r1).^2+2*A1*A2.*r1.*r2.*cos(k*(r1-r2)));
r1py = (Y-x1(2))./r1;                       r2py = (Y-x2(2))./r2;
ky = (k*(r1py.*A1^2.*r2.^2+ r2py.*A2^2.*r1.^2 + A1*A2*r1.*r2.*(r1py+r2py).*cos(k*(r1-r2)))...-
    - A1*A2*(r1py.*r2-r2py.*r1).*sin(k*(r1-r2)))./((A1*r2).^2+(A2*r1).^2+2*A1*A2.*r1.*r2.*cos(k*(r1-r2)));
K = sqrt(kx.^2+ky.^2)/k; 
dir = atan(kx./ky);

field = 1/(4*pi)*(A1*exp(-1i*k*r1)./r1 + A2*exp(-1i*k*r2)./r2);
    
A = sqrt((A1^2*r2.^2+A2^2*r1.^2+2*A1*A2*r1.*r2.*cos(k*(r1-r2)))./(r1.*r2));
%%
ftsize = 9;
fig = figure('Units','points','Position',[200,200,277,228]);
p1 = axes('Units','normalized','Position',[0.09 0.08 .9 .9 ]);

set(fig,'defaulttextinterpreter','latex')

pcolor(x/lambda,y/lambda,(K));
shading interp
axis equal tight
caxis([0 2]);
hold on
pcolor_ax =(gca);
cRange = caxis; 
[C,~] = contour( x/lambda, y/lambda , real(field), '-k');
hLines = findobj(gca, 'type', 'line');
set(hLines, 'LineWidth', 1);
caxis(cRange); 
xlabel( '$x/\lambda$', 'FontSize', ftsize );
ylabel( '$y/\lambda$', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
plot(y(1:5:end)*0,y(1:5:end)/lambda,'--k');
draw_ssd( fig, [x1;x2]/lambda, -[cosd(120) sind(120);cosd(60) sind(60)], 6e-2/lambda );
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);

c = colorbar;
title(c,'[dB]' , 'FontSize', ftsize-1);
set(gca,'FontName','Times New Roman');

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Appendices','stereo_wn.png' ) ,'-dpng')