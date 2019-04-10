clear
close all

fs =  8e3;
Lt =  0.65;
t = (-Lt:1/fs:Lt);
Nt = length(t);
w = 2*pi*(0:Nt/2)'/(Nt)*fs ;
c = 343.1;

dx = 0.03333;
dx_res = 0.10;
x0 = (-100:dx:100);
y0 = 0;
Nx = length(x0);
kx = 2*pi*(-Nx/2:Nx/2-1)'/(Nx)/dx;

xs = [ 0 -2 ];
yref = 1.5;

[X0,W] = meshgrid(x0,w);
r = sqrt( (X0-xs(1)).^2 + (y0-xs(2)).^2 );

D = -sqrt(1i*W/(2*pi*c)).*sqrt(yref/(yref-xs(2))).*xs(2).*exp(-1i*W/c.*r)./r.^(3/2);
Ds = zeros(size(D));
N0 = round(dx_res/dx);
Ds(1:N0:end) = D(1:N0:end)*N0;
Dkx = fftshift(fft(D,[],2),2)*dx;
Dkxs = fftshift(fft(Ds,[],2),2)*dx;
Dkx(isnan(Dkx))=0;
Dkxs(isnan(Dkxs))=0;

R = sqrt( X0.^2 + (yref-y0)^2 );
G0 = 1/(4*pi)*exp(-1i*W/c.*R)./R;
Gkx = fftshift(fft(G0,[],2),2)*dx;

R2 = sqrt( X0.^2 + (yref-xs(2))^2 );
G02 = 1/(4*pi)*exp(-1i*W/c.*R2)./R2;
Gkx2 = fftshift(fft(G02,[],2),2)*dx;
%%
kxs = 2*pi/dx_res;
f_a = c*pi/dx_res/(2*pi);
q = 12;
ftsize = 9;
fig = figure('Units','points','Position',[200,200,407,156]);
set(fig,'defaulttextinterpreter','latex')
colormap(flipud(pink))

pos = [ 0.06   0.17  0.38   .75
        0.53    0.17  0.485 .75];

p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(kx(1:q:end),w(1:q:end)/(2*pi)/1e3,20*log10(abs(Dkx(1:q:end,1:q:end))))
set(gca,'FontName','Times New Roman');
xlabel('$k_x$ [rad/m]')
ylabel('$f$ [kHz]')
shading interp
caxis([-30,10])
%
p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(kx(1:q:end),w(1:q:end)/(2*pi)/1e3,20*log10(abs(Dkxs(1:q:end,1:q:end))))
hold on
line([0,w(end)/c] ,[0,w(end)/(2*pi)/1e3],'Color','black','LineStyle','--')
line([0,-w(end)/c],[0,w(end)/(2*pi)/1e3],'Color','black','LineStyle','--')
%line([kx(1),kx(end)],[f_a,f_a]/1e3,'Color','black','LineStyle','.-')
shading interp
caxis([-30,10])
axis tight
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

xlabel('$k_x$ [rad/m]')
ylabel('$f$ [kHz]')

col = colorbar;
title(col,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-1);
%text(kxs-3,0.2,'$k_{x,s}$','FontSize',ftsize)
set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/SFS_theory','Aliased_spectrum' ) ,'-dpng')

%%
q = 12;
ftsize = 9;
fig = figure('Units','points','Position',[200,200,407,301]);
set(fig,'defaulttextinterpreter','latex')
colormap(flipud(pink))

pos = [ 0.065   0.6   0.38  .389
        0.59    0.6   0.38  .389
        0.33    0.085  0.38  .389];
 
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(kx(1:q:end)/kxs,w(1:q:end)/(2*pi)/1e3,20*log10(abs(Dkxs(1:q:end,1:q:end))))
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
shading interp
caxis([-30,10])
xlabel('$k_x/k_{x,s}$')
ylabel('$f$ [kHz]')

p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(kx(1:q:end)/kxs,w(1:q:end)/(2*pi)/1e3,20*log10(abs(Gkx(1:q:end,1:q:end))))
shading interp
axis tight
caxis([-90,-10])
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
xlabel('$k_x/k_{x,s}$')
ylabel('$f$ [kHz]');

%
p3 = axes('Units','normalized','Position',pos(3,:));
pcolor(kx(1:q:end)/kxs,w(1:q:end)/(2*pi)/1e3,20*log10(abs(Gkx(1:q:end,1:q:end).*Dkxs(1:q:end,1:q:end))))
shading interp
caxis([-90,-10])
axis tight
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
xlabel('$k_x/k_{x,s}$')
ylabel('$f$ [kHz]');


set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/SFS_theory','Aliased_repr_field' ) ,'-dpng')