clear 
close all

dx = 0.01;
x0 = (-50:dx:50)';
y0 = 0;
xs = [ 1.5,-2];


x = (-3:dx:3)';
y = (-1:dx:3)';
[X,Y] = meshgrid(x,y);

t0 = 0;

Lt =  1;
fs = 50e3;
t = (-Lt:1/fs:Lt);
Nt = length(t);
w = fftshift( 2*pi*(-Nt/2:Nt/2-1)'/(Nt)*fs );

w0 = 2*pi*1e3;
c = 343.1;
v = c*0.75;

kx = (w-w0)/v;
ky = -1i*sqrt( kx.^2 - (w/c).^2 );
yref = 1.5;

num   = besselh( 0, 2, ky*abs(yref-xs(2)));
denum = besselh( 0 ,2, ky*abs( yref-y0  ));
D0 = 1/v*num./denum.*tukeywin(Nt,0.0);
D0(isnan(D0)) = 0;

win = tukeywin(length(x0),.0);
%%

field_synth = zeros(size(X));

wb = waitbar(0,'Calculating field of SSD elements');
for n = 1 : length(x0)
   
    waitbar(n/length(x0),wb);
    R = sqrt( (X-x0(n)).^2 + (Y-y0).^2 );
    Tau = R/c;
 %   Tau_s = round(Tau*fs);
    
    Dwx = D0.*exp( -1i*kx*(x0(n)-xs(1)) );    
    dtx = win(n)*fftshift( ( ifft(Dwx) ) )*(w(2)-w(1))*Nt/(2*pi);
    S = interp1( t,dtx,t0 - Tau(:), 'spline' );
    field_synth = field_synth + 1/(4*pi)*reshape(S,length(y),length(x))./R*dx;
    
end
close(wb);
%%
M = v/c;
Delta_0 = sqrt( (X-xs(1)-v*t0).^2+ (Y-xs(2)).^2*(1-M^2) );
Tau_0 = (M*(X-xs(1)-v*t0) + Delta_0)/(c*(1-M^2));
field_ref = 1/(4*pi)*exp(1i*w0*(t0-Tau_0))./Delta_0;
%%

ftsize = 14.3;

% f = figure;
% set(f,'defaulttextinterpreter','latex')
% pcolor(x,y,real(field_ref));
% shading interp
% axis equal tight
% caxis([-1,1]*5e-2)
% xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
% ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
% set(gca,'FontName','Times New Roman');
% allAxesInFigure = findall(f,'type','axes');
% b = get(gca,'XTickLabel');
% set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
% line([x(1);x(end)],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
% line([x(1);x(end)],[0;0]+yref, 'Color', 'white','LineStyle','--','LineWidth',1);

fig = figure('Units','points','Position',[200,200,650,230]);
set(fig,'defaulttextinterpreter','latex')

pos = [ 0.06    0.09 0.38 .9
        0.53    0.09  0.465 .9];


p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field_synth) );
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
line([x(1);x(end)],[0;0]+yref, 'Color', 'white','LineStyle',':','LineWidth',1);


p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x,y,20*log10( abs( field_ref - field_synth ) ));
shading interp
axis equal tight
caxis([-75,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
col = colorbar;
title(col,'[dB]', 'FontSize', ftsize);

line([x(1);x(end)],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
line([x(1);x(end)],[0;0]+yref, 'Color', 'white','LineStyle',':','LineWidth',1);


%
set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','Linear_SDM' ) ,'-dpng')