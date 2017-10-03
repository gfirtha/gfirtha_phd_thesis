clear all
close all

dx = 0.01;
over_sample_L = 50;
L_ssd = 50;
x0 = (-over_sample_L*L_ssd:dx:over_sample_L*L_ssd)';
y0 = 0;
xs = [ 0,-2];
yref = 1.5;

Nx = length(x0);
kx = 2*pi*(-Nx/2:Nx/2-1)'/(Nx*dx);

f = 1000;
speed_of_sound = 343.1;
omega = 2*pi*f;
k = omega/speed_of_sound;
dkx = kx(2)-kx(1);
ky = -1i*sqrt(kx.^2 - k.^2);

D_kx = besselh( 0, 2, ky*abs(yref-xs(2)))./besselh( 0, 2, ky*abs(yref-y0))...
        .*exp(1i*kx*xs(1));

D_x_sdm = 1/(2*pi)*fftshift( fft( ifftshift( D_kx ) ) )*dkx.*tukeywin(length(x0),.1);
x_field = (-L_ssd:dx:L_ssd)';
D_x_sdm = spline(x0,D_x_sdm,x_field);
%
y_field = (-1:dx:3)';
[Xf,Yf] = meshgrid(x_field,y_field);
field = zeros(size(Xf));
%%
tic
wb = waitbar(0,'Calculating field of SSD elements');
for n = 1:length(y_field)
    
    waitbar(n/length(y_field),wb);
    r = sqrt( x_field.^2 + (y_field(n) - y0).^2 );
    G0 = 1/(4*pi)*exp(-1i*k*r)./r*dx;
    field(n,:) = conv(G0,D_x_sdm,'same');
    
    set( get(findobj(wb,'type','axes'),'title'), 'string', ...
                sprintf( 'Estimated time left %d seconds', round((length(y_field)-n)*toc/n)));
end
close(wb);
%%
r_ref = sqrt( (Xf-xs(1)).^2 + (Yf - xs(2)).^2 );
field_ref = 1/(4*pi)*exp(-1i*k*r_ref)./r_ref;
%%
L0 = 3;
x_fig = abs(x_field)<= L0;
ftsize = 13;
fig = figure('Units','points','Position',[200,200,650,230]);
pos = [ 0.06  0.065  0.38 .9
        0.53    0.065  0.465 .9];


p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x_field(x_fig),y_field,real(field(:,x_fig)));
shading interp
axis equal tight
caxis([-1,1]*5e-2)
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);
line([-L0;L0],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
line([-L0;L0],[0;0]+yref, 'Color', 'white','LineStyle',':','LineWidth',1);


p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x_field(x_fig),y_field,20*log10(abs(field_ref(:,x_fig)-field(:,x_fig))));
shading interp
axis equal tight
caxis([-85,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);
c = colorbar;
title(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);

line([-L0;L0],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
line([-L0;L0],[0;0]+yref, 'Color', 'white','LineStyle',':','LineWidth',1);


set(gcf,'PaperPositionMode','auto');
print( fullfile( '../..','Figures/SFS_theory','Linear_SDM' ) ,'-dpng')