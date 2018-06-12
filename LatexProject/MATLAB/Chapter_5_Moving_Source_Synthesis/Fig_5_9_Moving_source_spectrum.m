clear
close all


c = 343.1;
f0 = 0.25e3;
w0 = 2*pi*f0;
v = c*0.5;
M = v/c;
y = 1;
k = w0/c;

z = 0;
x = 0;
w = 2*pi*linspace(1,2*w0/(2*pi),100e3);

xs = 0;
ys = -2;
zs = 0;

t = linspace(-4e-2,4e-2,1e5 );
delta = sqrt( (x - xs - v*t).^2 + ( (y-ys)^2 + (z-zs)^2 )*(1-M^2) );
tau = ( M*(x-xs-v*t) +delta)/(c*(1-M^2));
pt = 1/(4*pi)*exp(1i*w0*(t-tau))./delta;

kx = ((w-w0)/v);
Pm = 1/v*(-1i/2*besselh(0,2, -1i*sqrt( kx.^2 - k.^2 ).*sqrt((y-ys)^2+(z-zs).^2))).*exp(-1i*kx.*(x-xs));

%%
ftsize = 14.3;
f = figure('Units','points','Position',[200,120,660,250]);
set(f,'defaulttextinterpreter','latex')

pos = [ 0.085  0.2 0.4 .72
        0.585 0.2 0.4 .72  ];


p1 = axes('Units','normalized','Position',pos(1,:));
plot(t*f0,real(pt),'LineWidth',1);
xlabel( '$t\cdot f_0 \rightarrow []$', 'FontSize', ftsize );
ylabel( '$P_{\mathrm{m}}(\mathbf{0},t,\omega_0) \rightarrow$ [Pa]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
grid on 

p2 = axes('Units','normalized','Position',pos(2,:));
plot(w/(2*pi)/f0 ,20*log10(abs(Pm)),'LineWidth',1);
ylim([-100,-30])
xlabel( '$f/f_0 \rightarrow []$', 'FontSize', ftsize );
ylabel( '$20\log_{10}|P_{\mathrm{m}}(\mathbf{0},\omega,\omega_0)| \rightarrow$ [dB]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
grid on

allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');

print( '-r300',fullfile( '../..','Figures/Moving_sources','moving_source_spectrum' ) ,'-dpng')