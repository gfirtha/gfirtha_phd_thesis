clear all
close all

dx = 2.5e-2;
L_ssd = 20;

x0_ = (-L_ssd :dx: L_ssd )';
y0 = 0;
z0_ =  (-L_ssd/2:dx: L_ssd/2 )';
[X0,Z0] = meshgrid(x0_, z0_);
Win = tukeywin(length(z0_),0.1)*tukeywin(length(x0_),0.1)';
Win = Win(:);
L_field = 3;
x_field = (-L_field:dx:L_field)';
y_field = (-1:dx:L_field)';
[ X_f, Y_f ] = meshgrid(x_field,y_field);
x0 = X0(:);
z0 = Z0(:);

xs = 1.5;
ys = -2;
zs = 0;
c = 343.1;
omega = 2*pi*1000;
k = omega/c;
v = c*0.75;
t = 0;
M = v / c;

%%
field_synth = zeros(size(X_f));
wb = waitbar(0,'Calculating field of SSD elements');
tic;


for n = 1:length(x0)
    waitbar(n/length(x0),wb);
    
    
    R0 = sqrt( ( X_f-x0(n)  ).^2 + ( Y_f ).^2 + z0(n).^2 );
    t_ret = t - R0/c;

    Delta = sqrt( (x0(n) - xs - v*t_ret).^2 + ((y0-ys).^2 + (z0(n)-zs).^2)*(1-M^2) ); 
    R = ( M * ( x0(n)-xs - v*t_ret ) + Delta ) /(1-M^2);
    Tau = R/c;
    Dx0 = 2*(y0-ys)/(4*pi)*( (1-M^2)./Delta + 1i*k )./ Delta.^2.*Win(n);

    
    field_synth = field_synth + 1/(4*pi)*Dx0.*exp( 1i * omega * (t - R0/c - Tau ))./R0*dx^2;

    
    set( get(findobj(wb,'type','axes'),'title'), 'string', ...
                sprintf( 'Estimated time left %d seconds', round((length(x0)-n)*toc/n)));
end
close(wb);
%%
Delta_0 = sqrt( (X_f-xs-v*t).^2+ (Y_f-ys).^2*(1-M^2) );
Tau_0 = (M*(X_f-xs-v*t) + Delta_0)/(c*(1-M^2));
field_ref = 1/(4*pi)*exp(1i*omega*(t-Tau_0))./Delta_0;
%%

ftsize = 9;
fig = figure('Units','points','Position',[200,200,407,144]);
set(fig,'defaulttextinterpreter','latex')

pos = [ 0.06    0.12  0.375   .84
        0.53    0.12  0.485  .84];



p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x_field,y_field,real(field_synth));
shading interp
axis equal tight
caxis([-1,1]*5e-2)
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
line([-L_field;L_field],[0;0], 'Color', 'black','LineStyle','-','LineWidth',2);

p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x_field,y_field,20*log10( abs( -field_ref + field_synth ) ));
shading interp
axis equal tight
caxis([-70,-10])
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
col = colorbar;
title(col,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-1);
line([-L_field;L_field],[0;0], 'Color', 'black','LineStyle','-','LineWidth',2);

set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/Moving_sources','Planar_SDM' ) ,'-dpng')