clear all
close all

R0 = 1.5;
fi = linspace(0,2*pi,180*2)';
dfi = fi(2)-fi(1);
x_ssd = R0*cos(fi);
y_ssd = R0*sin(fi);
n_ssd = [ -cos(fi) -sin(fi) ];

x = linspace(-2,2,300);
y = linspace(-2,2,300);
[X,Y] = meshgrid(x,y);

k = 2e3*2*pi/343.1;

% Generalized 2-1/2D WFS

x_s = -2.0;
y_s =  0;
v = [(x_ssd - x_s) (y_ssd - y_s)]; 
rx = v(:,1);
ry = v(:,2);
r0 = sqrt(dot(v,v,2));

Rref = 1.0;
Dr = R0^2 - Rref^2;
w = 2*(rx.*x_ssd + ry.*y_ssd);

A = Dr + r0.^2 - w;
B = r0.*(w - 2*Dr);
C = Dr.*r0.^2;

d_circ = ( -B - sqrt( B.^2 - 4*A.*C ) )./(2*A);
%%
w_ps = double ( dot(v,n_ssd,2)>0 );
D_ps = -sqrt(1i*k/(2*pi)).*(cos(fi).*v(:,1)+sin(fi).*v(:,2)).*...
        sqrt(d_circ).*exp(-1i*k*r0)./r0.^2.*w_ps;
    

%%
field_ps   = zeros(size(X));
for n = 1:length(x_ssd)
    r_ssd = sqrt(( X - x_ssd(n) ).^2 + ( Y - y_ssd(n) ).^2 );
    
    field_ps    = field_ps  + R0/(4*pi)*D_ps(n)...
                                            *exp(-1i*k*r_ssd)./r_ssd*dfi;
    
end

R = sqrt( (X-x_s).^2 + (Y-y_s).^2 );
field_ref_ps = 1/(4*pi)*exp(-1i*k*R)./R;

%%

ftsize = 15;
%%
g = figure('Units','points','Position',[50,50,750,350]);

subplot(1,2,1)
p1 = pcolor(x,y,20*log10(abs(field_ref_ps-field_ps)));

shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
caxis([-70,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,ps}}(\mathbf{x},\omega)-P_{\mathrm{ref,ps}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(g,'type','axes');
%b = get(gca,'XTickLabel');
%set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);
c = colorbar;
ylabel(c,'[dB]')
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);


subplot(1,2,2)
p2 = pcolor(x,y,real(field_ps));

shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
caxis([-.05,.05])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,ps}}(\mathbf{x},\omega)-P_{\mathrm{ref,ps}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(g,'type','axes');
%b = get(gca,'XTickLabel');
%set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);
c = colorbar;
ylabel(c,'[dB]')
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

set(gcf,'PaperPositionMode','auto');
