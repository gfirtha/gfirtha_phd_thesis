clear
close all

c = 343.1;
w = 2*pi*0.8e3;
k = w/c;

L = 2;
dx = 0.01;
x = (-L:dx:L)';
y = (-L:dx:L)';
[X,Y] = meshgrid(x,y);

r0 = sqrt(X.^2 + Y.^2);
p_ps = -1i/4*besselh(0,2,k*r0);

dfi = 0.1;
fi = (0:dfi:2*pi); R0 = 1.55;
x0_ps = R0*cos(fi); y0_ps = R0*sin(fi);
kx_ps = k*x0_ps./R0;
ky_ps = k*y0_ps./R0;
yref = 0.5;
nkx_ps2 = x./sqrt(x.^2+yref^2);

theta = 60;
x0_pw =  (-2*L:0.2:L)'; 
y0_pw = -x0_pw/tan(theta*pi/180);
kx_pw = k*cos(theta*pi/180)*ones(size(x0_pw));
ky_pw = k*sin(theta*pi/180)*ones(size(y0_pw));

p_pw = exp(-1i*k*(cos(theta*pi/180)*X+sin(theta*pi/180)*Y));
%%
ftsize = 12.2;
f = figure('Units','points','Position',[200,200,500,360]);

set(f,'defaulttextinterpreter','latex')
   
p1 = axes('Units','normalized','Position',[ 0.08 0.42 0.37 .6 ]);
pcolor(x,y,real(p_ps));
shading interp
axis equal tight
caxis([-.1,.1]);
hold on
contour( x, y, real(p_ps),[0 0], '-k');
plot(x0_ps,y0_ps,'--k')

headWidth = 5;
headLength = 5;
LineLength = 0.025;
for i = 1:1:length(x0_ps)
        ah = annotation('arrow',...
            'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[x0_ps(i) y0_ps(i) LineLength*kx_ps(i) LineLength*ky_ps(i)]);
end
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
xlabel( '$x \rightarrow$ [m]', 'FontSize', ftsize );
ylabel( '$y \rightarrow$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');

p2 = axes('Units','normalized','Position',[ 0.61 0.42 0.37 .6 ]);
pcolor(x,y,real(p_pw));
shading interp
axis equal tight
caxis([-1.5,1.5])
hold on
pcolor_ax =(gca);
cRange = caxis;
[C,~]= contour( x, y, real(p_pw),[0 0], '-k');
q = 50;

headWidth = 5;
headLength = 5;
LineLength = 0.025;
plot(x0_pw,y0_pw,'--k')
for i = 1:length(x0_pw)
        ah = annotation('arrow',...
            'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[x0_pw(i) y0_pw(i) LineLength*kx_pw(i) LineLength*ky_pw(i)]);
end
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

xlabel( '$x \rightarrow$ [m]', 'FontSize', ftsize );
ylabel( '$y \rightarrow$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
%
p3 = axes('Units','normalized','Position',[ 0.08 0.1 0.37 .2 ]);
plot(x,nkx_ps2);
xlabel( '$x \rightarrow$ [m]'  , 'FontSize', ftsize );
ylabel( '$\hat{k_x}^P(x,y_0) \rightarrow$'   , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
grid on

p4 = axes('Units','normalized','Position',[ 0.61 0.1 0.37 .2 ]);
plot(x,cos(theta*pi/180)*ones(size(x)));
xlabel( '$x \rightarrow$ [m]'  , 'FontSize', ftsize );
ylabel( '$\hat{k_x}^P(x,y_0) \rightarrow$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
ylim([0,1])
grid on

allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);


set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/High_freq_approximations','wavenumber_vector' ) ,'-dpng')