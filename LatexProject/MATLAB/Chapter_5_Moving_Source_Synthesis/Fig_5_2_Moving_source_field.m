clear
close all

c = 343.1;
w = 2*pi*1e3;
k = w/c;

dx = 1e-2;
x = (-2:dx:2)';
y = ( -2:dx:2)';

[X,Y] = meshgrid(x,y);
Z = 0;

t = 0.00;
v = c*0.5;
M = v/c;

xs = 0.0;
ys_e = 0;
zs_e = 0;

Delta = sqrt( ( X  - xs -v*t ).^2 + ((Y-ys_e).^2 + (Z-zs_e).^2)*(1-M^2) );
R = (M*( X-xs-v*t )+Delta)/(1-M^2);
Tau = R/c;

p = 1/(4*pi)*exp(1i*w*(t-Tau))./Delta;

xs_e0 = 0;
R0  = (xs-xs_e0)/M;
fi = linspace(0,2*pi,70)';
xr = [xs_e0+R0*cos(fi), 0+R0*sin(fi)];

kx = k*(xr(:,1)-xs_e0)./sqrt( ( xr(:,1) - xs -v*t ).^2 + xr(:,2).^2*(1-M^2) );
ky = k*(xr(:,2)-ys_e )./sqrt( ( xr(:,1) - xs -v*t ).^2 + xr(:,2).^2*(1-M^2) );
%%
ftsize = 9;
f = figure('Units','points','Position',[200,120,407,188]);
set(f,'defaulttextinterpreter','latex')

pos = [ 0.045 0.16 0.4 .76
        0.53  0.16 .48 .76  ];


p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(p));
shading interp
axis equal tight
hold on
contour( x, y, real(p),[0 0],'LineWidth',0.5,'LineColor',[1 1 1]*0.4);
caxis([-1,1]*.15);
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

headWidth = 5;
headLength = 5;
LineLength = 0.025;
plot(xr(:,1),xr(:,2),'--k')
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[xs 0 0.75 0]);
plot(xs,0,'ok','MarkerSize',2,'MarkerFaceColor','black')

p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x,y,20*log10(abs(p)));
shading interp
axis equal tight
caxis([-30,-20]);
hold on
contour( x, y, real(p),[0 0],'LineWidth',0.5,'LineColor',[1 1 1]*0.4);
%contour( x, y, abs(p),500,'LineWidth',0.5,'LineColor',[1 1 1]*0);


xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');

xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
c = colorbar;
title(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-1);

allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);


set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','moving_source_field' ) ,'-dpng')