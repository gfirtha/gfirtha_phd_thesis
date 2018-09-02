clear all
close

f = 1e3;
c = 343.1;
omega = 2*pi*f;
k = omega/c;

dx = 2.5e-2;
x0 = (-30:dx:30)';

x = linspace(-3,3,200);
y = linspace(-1,3,200);
[X,Y] = meshgrid(x,y);

xs = [1.5,-2];


yref = 1.5;

v = c*0.5;
t = 0;
M = v/c;
Dx0 = -sqrt(1i*k/(2*pi))*sqrt(yref/(yref-xs(2)))*xs(2);

field_synth = zeros( size( X ) );
wb = waitbar(0,'Calculating radiated field');

for n = 1 : length(x0)
    waitbar(n/length(x0),wb);
   
    R = sqrt( (X-x0(n)).^2 + (Y).^2 );
    t_ret = t - R/c;
    Delta = sqrt( (x0(n)-xs(1)-v*t_ret).^2 + (0-    xs(2)).^2*(1-M^2) );
    Tau = ( M * ( x0(n)-xs(1) - v*t_ret ) + Delta ) / ( c *(1-M^2) );

    field_synth = field_synth + 1/(4*pi)* Dx0*exp( 1i * omega * (t - R/c - Tau ))./(R.*Delta.^(3/2))*dx;
    
end
close(wb);
%%
Delta_0 = sqrt( (X-xs(1)-v*t).^2+ (Y-xs(2)).^2*(1-M^2) );
Tau_0 = (M*(X-xs(1)-v*t) + Delta_0)/(c*(1-M^2));
field_ref = 1/(4*pi)*exp(1i*omega*(t-Tau_0))./Delta_0;
%%
ftsize = 14.3;
fig = figure('Units','points','Position',[200,200,650,230]);
set(fig,'defaulttextinterpreter','latex')

pos = [ 0.06    0.09 0.38 .9
        0.53    0.09  0.465 .9];


p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field_synth));
shading interp
axis equal tight
caxis([-1,1]*5e-2)
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
line([x(1);x(end)],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
line([x(1);x(end)],[0;0]+yref, 'Color', 'white','LineStyle','--','LineWidth',1);


p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x,y,20*log10( abs( field_ref - field_synth ) ));
shading interp
axis equal tight
caxis([-73,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
col = colorbar;
title(col,'[dB]', 'FontSize', ftsize);

line([x(1);x(end)],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
line([x(1);x(end)],[0;0]+yref, 'Color', 'white','LineStyle','--','LineWidth',1);


%
set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/Moving_sources','25D_WFS_linear_SSD' ) ,'-dpng')