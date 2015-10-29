clear all
close all

c = 343.1;      % Speed of sound [m/s]
f = 2000;       % Plane wave frequency
w = 2*pi*f;     % Angular frequency
k = w/c;        % Acoustic wave number

x0 = [0 0 0];

dx = 0.005;
x = (-1:dx:1)';
y = (-1:dx:1)';
Z = 0;
[X,Y] = meshgrid(x,y);

r = sqrt( (X-x0(1)).^2 + (Y-x0(2)).^2 );
field_3D_greens = 1/(4*pi)*exp( 1i*k*r )./r;

cor = sqrt(1i*k/(2*pi));
field_2D_greens = -1i/4*besselh( 0, 2, k*r )*cor;
%%
ftsize = 8;
f = figure('Units','points','Position',[200,200,461,200]);
set(gcf,'Units','normalized');
subplot(1,2,1)
p1 = pcolor(x,y,real(field_3D_greens));
set(gca, 'Units','normalized','Position',[ 0.1 0.125 0.35 .9 ]);

shading interp
axis equal tight
caxis([-.25,.25]);
zlim([-1,1])
colormap gray

xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);

subplot(1,2,2)
p2 = plot( x, real( field_3D_greens(round(end/2),:) )  ,'k');
set(gca, 'Units','normalized','Position',[ 0.6 0.17 0.35 .812 ]);
ylim([-1.5,2])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$p \rightarrow [\mathrm{Pa} \hpace{2mm}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
%print -dpng point_source -r300
%% Dipole


r = sqrt( (X-x0(1)).^2 + (Y-x0(2)).^2 );
field_3D_dipole = -Y/(4*pi).*(1./r + 1i*w/c).*exp( 1i*k*r )./r.^2;

ftsize = 8;
f = figure('Units','points','Position',[200,200,230,200]);
set(gcf,'Units','normalized');
p1 = pcolor(x,y,real(field_3D_dipole));
set(gca, 'Units','normalized','Position',[ 0.1 0.125 .9 .85 ]);
shading interp
axis equal tight
caxis([-5,5]);
zlim([-1,1])
colormap gray

xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print -dpng dipole_source -r300