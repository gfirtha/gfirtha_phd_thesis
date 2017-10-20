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
ftsize = 11;
f = figure('Units','points','Position',[200,200,500,220]);
set(f,'defaulttextinterpreter','latex')

   
pos = [ 0.057  0.2 0.425 .8
        0.56  0.2 0.425 .8];
        
scale = 20;
q = 4;
p1 = axes('Units','normalized','Position',pos(1,:));
h0 = surf(y,y,scale*real(Gx));
%set(gca, 'Units','normalized','Position',[ 0.06 0.11 0.37 .9 ]);
shading interp
hold on
h1 =surf(x(1:q:end),y(1:q:end),scale*real(Gx(1:q:end,1:q:end)));
set(h1','edgecolor','k','facealpha',0,'linewidth',.1,'linestyle','-')
axis equal tight
caxis(scale*[-0.05,0.05]);
zlim([-2,3])
ax1 = gca;
ax1.ZTickLabel = ax1.ZTick/scale;
ax1.XTick = lambda*(-5:5:5);
ax1.YTick = lambda*(-5:5:5);

set(gca,'TickLabelInterpreter', 'tex');
tik_n = cell(5,1); tik_n{1} = '-5 \lambda';tik_n{2} = '0'; tik_n{3} = '5 \lambda';
ax1.XTickLabel = tik_n;
ax1.YTickLabel = tik_n;

xlabel('$x \rightarrow \mathrm{[m]}$','Interpreter','latex','FontSize',ftsize)
ylabel('$y \rightarrow \mathrm{[m]}$','Interpreter','latex','FontSize',ftsize)
zlabel('$\mathcal{R}\left(G(x,y,0,\omega_0)\right)$','Interpreter','latex','FontSize',ftsize)
set(gca,'FontName','Times New Roman');

q = 2;
scale = 10;
p2 = axes('Units','normalized','Position',pos(2,:));
h =surf(kx,ky,scale*abs(Gkx));
%set(gca, 'Units','normalized','Position',[ 0.528 0.11 0.47 .9 ]);
shading interp
hold on
h2 =surf(kx(1:q:end),ky(1:q:end),scale*abs(Gkx(1:q:end,1:q:end)));
set(h2','edgecolor','k','facealpha',0,'linewidth',.01,'linestyle',':')
axis equal
caxis([0,scale*10e-1]);
zlim([0,10])
ax2 = gca;
ax2.ZTickLabel = ax2.ZTick/scale;
ax2.XTick = k*(-1:1:1);
ax2.YTick = k*(-1:1:1);
tik_n2 = cell(7,1); 
tik_n2{1} = '-1\omega/\it{c}';
tik_n2{2} = '0';
tik_n2{3} = '1\omega/\it{c}';
ax2.XTickLabel = tik_n2;
ax2.YTickLabel = tik_n2;
set(gca,'TickLabelInterpreter', 'tex');
xlabel('$k_x \rightarrow \mathrm{[rad/m]}$','Interpreter','latex','FontSize',ftsize)
ylabel('$k_y \rightarrow \mathrm{[rad/m]}$','Interpreter','latex','FontSize',ftsize)
zlabel('$|(\tilde{G}(k_x,k_y,0,\omega_0)|$','Interpreter','latex','FontSize',ftsize);

set(gca,'FontName','Times New Roman');


set(gcf,'PaperPositionMode','auto');
print( '-r600', fullfile( '../..','Figures/Basic_acoustics','greens_function' ) ,'-dpng')