clear
close all

c = 343.1; 

fs = 3000;
T0 = 0.25;
t = (-T0:1/fs:T0)';
Nt = length(t);
w0 = 1000*2*pi;
w = 2*pi*(0:Nt-1)'/Nt*fs;

dx = 0.033;
Lx = 50;
x = (-Lx:dx:Lx)'; 
Nx = length(x);
kx = 2*pi*(-Nx/2:Nx/2-1)'/(Nx*dx);
[KX,W] = meshgrid(kx,w);

y = 0;
xs = 0; ys = -2;
yref = 1;

[X,T] = meshgrid(x,t);
v = c/2;
M = v/c;
input = zeros(size(t));

Wl = 1;
input(1 : Wl) = hann(Wl);
input_sp = (fft(input));

R_st = sqrt((X-xs).^2 + (y-ys).^2);

Delta = sqrt( ( (X-xs) - v*T ).^2 + (y-ys)^2*(1-M^2) );
R_dyn = ( M*( X-xs-v*T ) + Delta  )/(1-M^2);

P_st = 1/(4*pi)*exp( 1i*w0*( T - R_st/c ) )./R_st;
P_dy = 1/(4*pi)*exp( 1i*w0*( T - R_dyn/c) )./Delta;

%
Win = tukeywin(Nt,.5)*tukeywin(Nx,.5)';
Pkx_stationary = fftshift( fft2( P_st.*Win ) , 2 );
Pkx_moving     = fftshift( fft2( P_dy.*Win ) , 2 );
%%
q = 2;
ftsize = 9;
fig = figure('Units','points','Position',[200,200,407,156]);
set(fig,'defaulttextinterpreter','latex')
colormap(flipud(pink))

pos = [ 0.07   0.17  0.38   .75
        0.55    0.17  0.465 .75];
cax = [-50,50];
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(kx(1:q:end),w(1:q:end)/w0,20*log10(abs(Pkx_stationary(1:q:end,1:q:end))))
hold on
line([0,w(end)/c],[0,w(end)/w0],'Color','black','LineStyle','--')
line([0,-w(end)/c],[0,w(end)/w0],'Color','black','LineStyle','--')
line([0,0],[0,w(end)]/w0,'Color','black','LineStyle',':')
set(gca,'FontName','Times New Roman');
% 

xlabel('$k_x$ [rad/m]')
ylabel('$f/f_0$')
xlim([-70,70])
shading interp
caxis(cax)
%
p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(kx(1:q:end),w(1:q:end)/w0,20*log10(abs(Pkx_moving(1:q:end,1:q:end))))
hold on
line([0,w(end)/c],[0,w(end)/w0],'Color','black','LineStyle','--')
line([0,-w(end)/c],[0,w(end)/w0],'Color','black','LineStyle','--')
line([0,0],[0,w(end)]/w0,'Color','black','LineStyle',':')
line([kx(1),kx(end)],[w0 w0]/w0,'Color','black','LineStyle',':')
shading interp
caxis(cax)
axis tight
set(gca,'FontName','Times New Roman');
xlabel('$k_x$ [rad/m]')
ylabel('$f/f_0$')
xlim([-70,70])
col = colorbar;
title(col,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-1);

allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','Moving_source_kxw' ) ,'-dpng')