clear all
close

f = 1e3;
speed_of_sound = 343.1;
omega = 2*pi*f;
k = omega/speed_of_sound;

dx = 1e-2;
x0 = (-20:dx:20)';

x = linspace(-3,3,200);
y = linspace(-1,3,200);
[X,Y] = meshgrid(x,y);

xs = [0,-2];

r0 = sqrt( (x0-xs(1)).^2 + (-xs(2)).^2 );
yref = 1.5;

Dx0 = -1/(4*pi)*sqrt( 8*pi/(1i*k) )*sqrt(yref/(yref-xs(2)))*1i*k*xs(2)*exp(-1i*k*r0)./r0.^(3/2);

field_synth = zeros( size( X ) );
wb = waitbar(0,'Calculating radiated field');

for n = 1 : length(x0)
    waitbar(n/length(x0),wb);
   
    R = sqrt( (X-x0(n)).^2 + (Y).^2 );
    field_synth = field_synth + 1/(4*pi)*Dx0(n)*exp(-1i*k*R)./R*dx;
    
end
close(wb);

field_ref = 1/(4*pi)*exp(-1i*k*sqrt( (X-xs(1)).^2 + (Y-xs(2)).^2 ))./sqrt( (X-xs(1)).^2 + (Y-xs(2)).^2 );
%%
ftsize = 13;
fig = figure('Units','points','Position',[200,200,650,230]);
pos = [ 0.06  0.065  0.38 .9
        0.53    0.065  0.465 .9];


p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field_synth));
shading interp
axis equal tight
caxis([-1,1]*5e-2)
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);
line([x(1);x(end)],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
line([x(1);x(end)],[0;0]+yref, 'Color', 'white','LineStyle',':','LineWidth',1);


p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x,y,20*log10( abs( field_ref - field_synth ) ));
shading interp
axis equal tight
caxis([-73,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);
c = colorbar;
title(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);

line([x(1);x(end)],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
line([x(1);x(end)],[0;0]+yref, 'Color', 'white','LineStyle',':','LineWidth',1);


%
set(gcf,'PaperPositionMode','auto');
print( fullfile( '../..','Figures/SFS_theory','25D_WFS_linear_SSD' ) ,'-dpng')