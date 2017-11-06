clear 
close all

c = 343.1;
w = 2*pi*1e3;
k = w/c;

dx = 0.025;
x0 = (-10:dx:10)';
y0 = 0;

x_field = (-3:dx:3)';
y_field = (-1:dx:3)';
[X,Y] = meshgrid( x_field, y_field );
Rref = 4;
ys = -2;
xs = 0;

r0 = sqrt((x0-xs).^2+(y0-ys)^2);
xref = zeros(length(x0),2);
xref(:,1) = xs + (x0-xs).*Rref./r0;
xref(:,2) = ys + (y0-ys).*Rref./r0;
%
R = sqrt( (xref(:,1)-xs).^2 + (xref(:,2) - ys).^2 );
dref = Rref-r0;

G = 1/(4*pi)*exp( -1i*k*dref )./dref;

P_2D = -1i/4*besselh(0,2, k*R);
P_3D = 1/(4*pi)*exp( -1i*k*R )./R;

%ky = (xref(:,2)-ys)./R;
%Dx_3D =sqrt( R./(R-dref) ).*sqrt(1i*k./(2*pi*dref)).*ky.*P_3D./G;
Dx_3D = (-ys)*sqrt(1i*k/(2*pi))*sqrt((Rref-r0)./Rref).*exp(-1i*k*r0)./r0.^(3/2);

Dx_3D(abs(x0-xs)>=sqrt(Rref^2-ys^2)) = 0;
win = zeros(length(Dx_3D),1);
win(Dx_3D~=0) = tukeywin(sum(double(Dx_3D~=0)),0.2);
Dx_3D = Dx_3D.*win;
%%
field_synth_3D = zeros(size(X));

wb = waitbar(0,'Calculating ');
for n = 1:length(x0)
    waitbar(n/length(x0),wb);
    R0 = sqrt( (X-x0(n)).^2 + Y.^2 );
    G0_3D = 1/(4*pi)*exp(-1i*k*R0)./R0;
    field_synth_3D = field_synth_3D + Dx_3D(n)*G0_3D*dx;
end
close(wb)

R = sqrt((X-xs).^2 + (Y-ys).^2);
field_ref_3D = 1/(4*pi)*exp(-1i*k*R)./R;
%%

ftsize = 14.3;
fig = figure('Units','points','Position',[200,200,650,230]);
set(fig,'defaulttextinterpreter','latex')

pos = [ 0.06    0.09 0.38 .9
        0.53    0.09  0.465 .9];

p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x_field,y_field,real(field_synth_3D));
shading interp
axis equal tight
%caxis([-1,1]*5e-2)

caxis([-.05,.05])
hold on
plot( [ -3+dx 3-dx ],  [ 0 0 ], 'k', 'Linewidth', 2 )
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);
line([x_field(1);x_field(end)],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);

lim = x0(abs(x0)<sqrt(Rref^2-ys^2));
xc = lim;
yc = sqrt(Rref^2 - xc.^2) + ys;

p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x_field,y_field,20*log10(abs(field_ref_3D-field_synth_3D)));
shading interp
axis equal tight
hold on
plot( [ -3+dx 3-dx ],  [ 0 0 ], 'k', 'Linewidth', 2 )
caxis([-70,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
c = colorbar;
title(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);
plot(xc,yc, 'Color', 'white','LineStyle',':','LineWidth',1);
line([x_field(1);x_field(end)],[0;0], 'Color', 'black','LineStyle','-','LineWidth',1);
xlim([x_field(1),x_field(end)]);

%
set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/SFS_theory','25D_spatial_SDM_linear_SSD' ) ,'-dpng')