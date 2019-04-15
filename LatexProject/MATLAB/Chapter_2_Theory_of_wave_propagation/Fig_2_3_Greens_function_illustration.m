clear
close all

c = 343;
w = 2*pi*500;
k = w/c * (1-0.00*1i*1e-2);
x0  = 0;   y0 = 0;    z0 = 0;

lambda= 2*pi*c/w;

dx = 0.025;
L = 5*lambda;
x = (-L:dx:L+dx)';  y = x;
Nx = length(x);

kxs = 1.5*k;
dkx = kxs/500;
kx =   (-kxs:dkx:kxs )';
ky =  kx;
dkx = kx(2)-kx(1);
[Kx,Ky] = meshgrid(kx,ky);
Gkx = -1i/2*exp(1i*(Kx*x0+Ky*y0)).*exp(-sqrt(-k^2 + (Kx.^2 + Ky.^2))*abs(z0))./ sqrt(k^2 - (Kx.^2 + Ky.^2));

[X,Y] = meshgrid(x,y);
Gx = 1/(4*pi)*exp(-1i*k*sqrt((X-x0).^2+(Y-y0).^2 + z0.^2))./sqrt((X-x0).^2+(Y-y0).^2 + z0.^2);
%%
ftsize = 9;
f = figure('Units','points','Position',[200,200,407,170]);
set(f,'defaulttextinterpreter','latex')

   
pos = [ 0.035   0.2 0.45 .75
        0.575  0.2 0.45 .75];
        
scale = 20;
q = 4;
p1 = axes('Units','normalized','Position',pos(1,:));
h0 = surf(x/lambda,y/lambda,scale*real(Gx));
%set(gca, 'Units','normalized','Position',[ 0.06 0.11 0.37 .9 ]);
shading interp
hold on
h1 =surf(x(1:q:end)/lambda,y(1:q:end)/lambda,scale*real(Gx(1:q:end,1:q:end)));
set(h1','edgecolor','k','facealpha',0,'linewidth',.1,'linestyle','-')
axis equal tight
caxis(scale*[-0.05,0.05]);
zlim([-2,3])

set(gca,'TickLabelInterpreter', 'tex');

xlabel('$x/\lambda$','Interpreter','latex','FontSize',ftsize)
ylabel('$y/\lambda$','Interpreter','latex','FontSize',ftsize)
zlabel('$\mathcal{R}\left(G(x,y,0,\omega_0)\right)$','Interpreter','latex','FontSize',ftsize)
set(gca,'FontName','Times New Roman');

q = 2;
p2 = axes('Units','normalized','Position',pos(2,:));
h =surf(kx/k,ky/k,abs(Gkx));
shading interp
hold on
h2 =surf(kx(1:q:end)/k,ky(1:q:end)/k,abs(Gkx(1:q:end,1:q:end)));
set(h2','edgecolor','k','facealpha',0,'linewidth',.01,'linestyle',':')
axis equal
caxis([0,1]*10e-1);
zlim([0,1])
set(gca,'TickLabelInterpreter', 'tex');
xlabel('$k_x/k$','Interpreter','latex','FontSize',ftsize)
ylabel('$k_y/k$','Interpreter','latex','FontSize',ftsize)
zlabel('$|(\tilde{G}(k_x,k_y,0,\omega_0)|$','Interpreter','latex','FontSize',ftsize);

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);


set(gcf,'PaperPositionMode','auto');
print( '-r600', fullfile( '../..','Figures/Basic_acoustics','greens_function' ) ,'-dpng')