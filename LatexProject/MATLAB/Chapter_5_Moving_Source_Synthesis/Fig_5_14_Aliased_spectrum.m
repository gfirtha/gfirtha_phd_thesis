clear
close all

c = 343.1;

fs = 4000;
T0 = 0.1;
t = (-T0:1/fs:T0)';
Nt = length(t);
%w0 = 1000*2*pi;
w0 = linspace(0,2000,200)*2*pi;
Hw0 = tukeywin(2*length(w0),0.1);
Hw0 = Hw0(end/2+1:end);
w = 2*pi*(0:Nt-1)'/Nt*fs;

dx = 0.03333;
dx_res = 0.10;
Lx = 50;
x = (-Lx:dx:Lx)';
Nx = length(x);
kx = 2*pi*(-Nx/2:Nx/2-1)'/(Nx*dx);
[KX,W] = meshgrid(kx,w);

yref = 1.5;
xs = 0; ys = -2;

[X,T] = meshgrid(x,t);
v = c*2/4;
M = v/c;
input = zeros(size(t));

Wl = 1;
input(1 : Wl) = hann(Wl);
input_sp = (fft(input));


Delta = sqrt( ( (X-xs) - v*T ).^2 + (yref-ys)^2*(1-M^2) );
R_dyn = ( M*( X-xs-v*T ) + Delta  )/(1-M^2);
D = 0;
tic
for i = 1  : length(w0)
    i
D = D-Hw0(i)*sqrt(-1i*w0(i)/(2*pi*c))*sqrt(yref/(yref-ys))*ys*exp( 1i*w0(i)*( T - R_dyn/c) )./Delta.^(3/2)*mean(diff(w0));
end
toc
Ds = zeros(size(D));
N0 = round(dx_res/dx);
Ds(:,1:N0:end) = D(:,1:N0:end)*N0;
%
Win = tukeywin(Nt,.5)*tukeywin(Nx,.5)';
dw = (2*pi)/fs;
Dkx = fftshift( fft2( D.*Win ) , 2 )*dx/fs;
Dkxs = fftshift( fft2( Ds.*Win ) , 2 )*dx/fs;

R = sqrt( X.^2 + (yref)^2 );
G0 = 1/(4*pi)*exp(-1i*W/c.*R)./R;
Gkx = fftshift(fft(G0,[],2),2)*dx;
%%
kxs = 2*pi/dx_res;
q = 2;
ftsize = 14.3;
fig = figure('Units','points','Position',[200,200,730,540]);
set(fig,'defaulttextinterpreter','latex')
colormap(flipud(pink))

pos = [ 0.065   0.59   0.38  .389
        0.59    0.59   0.38  .389
        0.33    0.075  0.38  .389];
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(kx(1:q:end)/kxs,w(1:q/2:end)/(2*pi*1e3),20*log10(abs(Dkxs(1:q/2:end,1:q:end))))
hold on
line([0,w(end)/c/kxs],[0,w(end)/(2*pi*1e3)],'Color','black','LineStyle','--')
line([0,-w(end)/c/kxs],[0,w(end)/(2*pi*1e3)],'Color','black','LineStyle','--')
%line([0,0],[0,w(end)/(2*pi*1e3)],'Color','black','LineStyle',':')
set(gca,'FontName','Times New Roman');
xlabel('$k_x/k_{x,s} \rightarrow$ []')
ylabel('$f \rightarrow$ [kHz]');
shading interp
caxis([-50,30])
axis tight
%xlim([-70,70]/kxs)
%ylim([0,5]);

p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(kx(1:q:end)/kxs,w(1:q/2:end)/(2*pi)/1e3,20*log10(abs(Gkx(1:q/2:end,1:q:end))))
shading interp
axis tight
caxis([-90,-10])
hold on
line([0, w(end)/c/kxs],[0,w(end)/(2*pi*1e3)],'Color','black','LineStyle','--')
line([0,-w(end)/c/kxs],[0,w(end)/(2*pi*1e3)],'Color','black','LineStyle','--')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
xlabel('$k_x/k_{x,s} \rightarrow$ []')
ylabel('$f \rightarrow$ [kHz]');
%xlim([-70,70]/kxs)
%ylim([0,5]);

p3 = axes('Units','normalized','Position',pos(3,:));
pcolor(kx(1:q:end)/kxs,w(1:q/2:end)/(2*pi*1e3),20*log10(abs(Gkx(1:q/2:end,1:q:end).*Dkxs(1:q/2:end,1:q:end))))
hold on
line([0,w(end)/c/kxs],[0,w(end)/(2*pi*1e3)],'Color','black','LineStyle','--')
line([0,-w(end)/c/kxs],[0,w(end)/(2*pi*1e3)],'Color','black','LineStyle','--')
%line([0,0],[0,w(end)/(2*pi*1e3)],'Color','black','LineStyle',':')
shading interp
cax = [-80,-5];
caxis(cax)
axis tight
set(gca,'FontName','Times New Roman');
xlabel('$k_x/k_{x,s} \rightarrow$ []')
ylabel('$f \rightarrow$ [kHz]');
%xlim([-70,70]/kxs)
%ylim([0,5]);
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','Aliased_spectrum' ) ,'-dpng');