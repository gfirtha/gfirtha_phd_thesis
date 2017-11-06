clear all
close all

dx = 0.025;
L_ssd = 15;

x0 = (-L_ssd :dx: L_ssd )';
y0 = 0;
z0 = (-L_ssd :dx :L_ssd )';
[X0,Z0] = meshgrid(x0, z0);
Win = tukeywin(length(z0),0.1)*tukeywin(length(x0),0.1)';

L_field = 3;
x_field = (-L_ssd:dx:L_ssd)';
y_field = (-1:dx:L_field)';
[ X_f, Y_f ] = meshgrid(x_field,y_field);

xs = [0, -2, 0];
c = 343.1;
w = 2*pi*1000;
k = w/c;

R0 = sqrt( (X0 - xs(1)).^2 + (y0-xs(2)).^2 + (Z0-xs(3)).^2 ); 
Dx = 2*(y0-xs(2))/(4*pi)* ( 1./R0 + 1i*k ).*exp( -1i*k*R0 ) ./ R0.^2.*Win;

%%
field = zeros(size(X_f));

wb = waitbar(0,'Calculating field of SSD elements');
tic;
for z = 1:length(z0)
    waitbar(z/length(z0),wb);
    for y = 1:length(y_field) 
    R = sqrt( x_field.^2 + ( y_field(y)-y0 ).^2 + z0(z).^2 );
    G0 = 1/(4*pi)*exp(-1i*k*R)./R*dx^2;
    field(y,:) = field(y,:) + conv( G0, Dx(:,z),'same' ).';
    end
    set( get(findobj(wb,'type','axes'),'title'), 'string', ...
                sprintf( 'Estimated time left %d seconds', round((length(z0)-z)*toc/z)));
end
close(wb);
%%
r_ref = sqrt( (X_f-xs(1)).^2 + (Y_f-xs(2)).^2 + (xs(3)).^2 );
ref_field = 1/(4*pi)*exp(-1i*k*r_ref)./r_ref;
%%
x_fig = abs(x_field)<= 3;
dx = x_field(2)-x_field(1);

ftsize = 14.3;
fig = figure('Units','points','Position',[200,200,650,230]);
set(fig,'defaulttextinterpreter','latex')

pos = [ 0.06    0.09 0.38 .9
        0.53    0.09  0.465 .9];


p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x_field(x_fig),y_field,real(field(:,x_fig)));
shading interp
axis equal tight
caxis([-1,1]*3e-2)
xlabel( '$x \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
line([-L_field;L_field],[0;0], 'Color', 'black','LineStyle','-','LineWidth',2);


p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x_field(x_fig),y_field,20*log10( abs( ref_field(:,x_fig) - field(:,x_fig) ) ));
shading interp
axis equal tight
caxis([-70,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
c = colorbar;
title(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);

line([-L_field;L_field],[0;0], 'Color', 'black','LineStyle','-','LineWidth',2);



set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/SFS_theory','Planar_SDM' ) ,'-dpng')