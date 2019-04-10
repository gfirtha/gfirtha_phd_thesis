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
yref = 2;
[X0,W] = meshgrid(x0,w);
r = sqrt( (X0-xs(1)).^2 + (y0-xs(2)).^2 );

kxs = 2*pi/dx_res;
D = -sqrt(1i*W/(2*pi*c)).*sqrt(yref/(yref-xs(2))).*xs(2).*exp(-1i*W/c.*r)./r.^(3/2);
khx0 = (X0-xs(1))./r;
KxNy = pi/dx_res;
kx0 = -0.25;
A = get_transfer( W, KxNy*c./abs(khx0-kx0), .75);

D = A.*D;
%D( W >= (KxNy*c./abs(khx0)) ) = 0;
Ds = zeros(size(D));
N0 = round(dx_res/dx);
Ds(1:N0:end) = D(1:N0:end)*N0;
Dkx = fftshift(fft(D,[],2),2)*dx;
Dkxs = fftshift(fft(Ds,[],2),2)*dx;
Dkx(isnan(Dkx))=0;
Dkxs(isnan(Dkxs))=0;
%%
q = 10;
ftsize = 9;
fig = figure('Units','points','Position',[200,200,260,144]);
set(fig,'defaulttextinterpreter','latex')
colormap(flipud(pink))

pos = [ 0.11   0.22  0.8   .75];
%
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(kx(1:q:end)/kxs,w(1:q:end)/(2*pi)/1e3,20*log10(abs(Dkxs(1:q:end,1:q:end))))
hold on
line([-1-kx0 -1-kx0]/2,[0,w(end)/(2*pi)/1e3],'Color','black','LineStyle','--','LineWidth',1)
line([1-kx0 1-kx0]/2,[0,w(end)/(2*pi)/1e3],'Color','black','LineStyle','--','LineWidth',1)
line(-[kx0 kx0]/2,[0,w(end)/(2*pi)/1e3],'Color','black','LineStyle',':','LineWidth',1)

shading interp
caxis([-30,5])
axis tight
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
xlabel('$k_x/k_{x,s}$')
ylabel('$f$ [kHz]')

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/SFS_theory','AntiAliased_spectrum_asymm' ) ,'-dpng')
