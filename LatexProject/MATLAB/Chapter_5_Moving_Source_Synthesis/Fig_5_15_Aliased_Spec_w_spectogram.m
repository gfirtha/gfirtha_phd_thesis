clear
close all
addpath(genpath('../Files'));

c = 343.1;

fs = 3000;
T0 = 0.25;
t = (-T0:1/fs:T0)';
Nt = length(t);
w0 = 1000*2*pi;
w = 2*pi*(0:Nt-1)'/Nt*fs;

dx = 0.03333;
dx_res = 0.1;
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
D = -sqrt(-1i*w0/(2*pi*c))*sqrt(yref/(yref-ys))*ys*exp( 1i*w0*( T - R_dyn/c) )./Delta.^(3/2);
Ds = zeros(size(D));
N0 = round(dx_res/dx);
Ds(:,1:N0:end) = D(:,1:N0:end)*N0;
Win = tukeywin(Nt,.5)*tukeywin(Nx,.5)';
dw = (2*pi)/fs;
Dkx = fftshift( fft2( D.*Win ) , 2 )*dx/fs;
Dkxs = fftshift( fft2( Ds.*Win ) , 2 )*dx/fs;
R = sqrt( X.^2 + (yref)^2 );
G0 = 1/(4*pi)*exp(-1i*W/c.*R)./R;
Gkx = fftshift(fft(G0,[],2),2)*dx;

Pkx = Gkx.*Dkxs;
%%
Lx = 50;
x = (-Lx:dx_res:Lx)';
Pxt = 0;

fs = 44100;
T0 = 0.4;
t = (-T0:1/fs:T0)';
for n = 1 : length(x)
    r0 = sqrt( x(n).^2 + yref.^2 );
    Delta = sqrt( ( (x(n)-xs) - v*(t-r0/c) ).^2 + (yref-ys)^2*(1-M^2) );
    R_dyn = ( M*( x(n)-xs-v*(t-r0/c)  ) + Delta  )/(1-M^2);
    Pxt = Pxt + 1/(4*pi)*1./r0*dx_res.*...
        -sqrt(-1i*w0/(2*pi*c))*sqrt(yref/(yref-ys))*ys*exp( 1i*w0*( t - r0/c - R_dyn/c) )./Delta.^(3/2);
end

Win_length = 5e-2;
step = Win_length/30;
wu = w0/(1-M)*2;
wd = w0/(1+M)*0;
w_ = linspace( wd, wu, 5000 )';
f_= w_/(2*pi);
%
[ stft_synth, t_]     =  stft( Pxt.' , w_, fs, Win_length, step ) ;
%%

kxs = 2*pi/dx_res;
q = 2;
ftsize = 9;
fig = figure('Units','points','Position',[200,200,407,156]);
set(fig,'defaulttextinterpreter','latex')
colormap(flipud(pink))

pos = [ 0.08  0.17  0.38   .75
        0.57  0.17  0.38 .75];
cax = [-130,-30];
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(kx(1:q:end)/kxs,w(1:q:end)/w0,20*log10(abs(Pkx(1:q:end,1:q:end))))
hold on
line([0,w(end)/c/kxs],[0,w(end)/w0],'Color','black','LineStyle','--')
line([0,-w(end)/c/kxs],[0,w(end)/w0],'Color','black','LineStyle','--')
line([0,0],[0,w(end)]/w0,'Color','black','LineStyle',':')
line([kx(1),kx(end)]/kxs,[w0 w0]/w0,'Color','black','LineStyle',':')
set(gca,'FontName','Times New Roman');
xlabel('$k_x/k_{x,s}$')
ylabel('$f/f_0$')
shading interp
cax = [-130,-30]-20;
caxis(cax)
axis tight
xlim([-70,70]/kxs)
ylim([0,3])


p2 = axes('Units','normalized','Position',pos(2,:));
pcolor( t_+t(1)-step*7, f_/w0*2*pi, (abs(stft_synth)) );
hold on
% line([0,w(end)/c/kxs],[0,w(end)/w0],'Color','black','LineStyle','--')
% line([0,-w(end)/c/kxs],[0,w(end)/w0],'Color','black','LineStyle','--')
line([0,0],[0,w(end)]/w0,'Color','black','LineStyle',':')
line([(t_(1)-T0),(t_(end)-T0)],[w0 w0]/w0,'Color','black','LineStyle',':')
set(gca,'FontName','Times New Roman');
xlabel('$t$ [s]')
ylabel('$f/f_0$ ')
shading interp
%cax = [-110,-75];
caxis([0 1]*8e-5)
axis tight
ylim([0,3])
xlim([-0.25,0.25])


allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','Aliased_spectogram' ) ,'-dpng');