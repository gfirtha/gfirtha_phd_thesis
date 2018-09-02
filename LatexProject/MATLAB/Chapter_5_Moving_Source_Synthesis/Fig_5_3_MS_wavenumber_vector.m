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

xs = 0.4;
ys_e = 0;
zs_e = 0;

Delta = sqrt( ( X  - xs -v*t ).^2 + ((Y-ys_e).^2 + (Z-zs_e).^2)*(1-M^2) );
R = (M*( X-xs-v*t )+Delta)/(1-M^2);
Tau = R/c;

p = 1/(4*pi)*exp(1i*w*(t-Tau))./Delta;

xs_e0 = -0.155;
R0  = (xs-xs_e0)/M;
fi = linspace(0,2*pi,70)';
xr = [xs_e0+R0*cos(fi), 0+R0*sin(fi)];

kx = k*(xr(:,1)-xs_e0)./sqrt( ( xr(:,1) - xs -v*t ).^2 + xr(:,2).^2*(1-M^2) );
ky = k*(xr(:,2)-ys_e )./sqrt( ( xr(:,1) - xs -v*t ).^2 + xr(:,2).^2*(1-M^2) );
%%
ftsize = 13;
f = figure('Units','points','Position',[200,120,320,300]);
set(f,'defaulttextinterpreter','latex')

% pos = [ 0.1 0.16 0.86 .86 ];
% 
%    
% p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(p));
shading interp
axis equal tight

hold on
contour( x, y, real(p),[0 0],'-k');
caxis([-1,1]*.15);
xlabel( '$x \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
plot(xs_e0,0,'ok','MarkerSize',3,'MarkerFaceColor','black')

headWidth = 5;
headLength = 5;
LineLength = 0.025;
plot(xr(:,1),xr(:,2),'--k')
for i = 1:length(xr)
        ah = annotation('arrow',...
            'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[xr(i,1) xr(i,2) LineLength*kx(i) LineLength*ky(i)]);
end
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);


set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','moving_source_local_wave' ) ,'-dpng')