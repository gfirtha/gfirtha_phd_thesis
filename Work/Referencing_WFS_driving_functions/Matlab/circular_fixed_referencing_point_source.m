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
d_fixed = 0.75;

x_s = -2.0;
y_s =  0;
v = [(x_ssd - x_s) (y_ssd - y_s)]; 
r0 = sqrt(dot(v,v,2));

w_ps = double ( dot(v,n_ssd,2)>0 );
D_ps = -sqrt(1i*k/(2*pi)).*(cos(fi).*v(:,1)+sin(fi).*v(:,2)).*...
        sqrt(d_fixed).*exp(-1i*k*r0)./r0.^2.*w_ps;
    
    
x_s_2 = -3.0;
v_2 = [(x_ssd - x_s_2) (y_ssd - y_s)]; 
r0_2 = sqrt(dot(v_2,v_2,2));
w_ps_2 = double ( dot(v_2,n_ssd,2)>0 );
D_ps_2 = -sqrt(1i*k/(2*pi)).*(cos(fi).*v_2(:,1)+sin(fi).*v_2(:,2)).*...
        sqrt(d_fixed).*exp(-1i*k*r0_2)./r0_2.^2.*w_ps_2;
%%
field_ps   = zeros(size(X));
field_ps_2  = zeros(size(X));
for n = 1:length(x_ssd)
    r_ssd = sqrt(( X - x_ssd(n) ).^2 + ( Y - y_ssd(n) ).^2 );
    
    field_ps    = field_ps  + R0/(4*pi)*D_ps(n)...
                                            *exp(-1i*k*r_ssd)./r_ssd*dfi;
    field_ps_2  = field_ps_2  + R0/(4*pi)*D_ps_2(n)...
                                            *exp(-1i*k*r_ssd)./r_ssd*dfi;
    
end

R = sqrt( (X-x_s).^2 + (Y-y_s).^2 );
field_ref_ps = 1/(4*pi)*exp(-1i*k*R)./R;

R2 = sqrt( (X-x_s_2).^2 + (Y-y_s).^2 );
field_ref_ps_2 = 1/(4*pi)*exp(-1i*k*R2)./R2;
%%
%
v_n = bsxfun( @times, v, 1./r0);
alpha_0 = acos(v_n(:,1)).*sign(v_n(:,2));
%
x_ps = x_s + (r0./(r0-d_fixed)).*cos(alpha_0).*r0;
y_ps = y_s + (r0./(r0-d_fixed)).*sin(alpha_0).*r0;

w_ps(w_ps==0) = nan;

v_n2 = bsxfun( @times, v_2, 1./r0_2);
alpha_0_2 = acos(v_n2(:,1)).*sign(v_n2(:,2));
%
x_ps2 = x_s_2 + (r0_2./(r0_2-d_fixed)).*cos(alpha_0_2).*r0_2;
y_ps2 = y_s   + (r0_2./(r0_2-d_fixed)).*sin(alpha_0_2).*r0_2;

w_ps_2(w_ps_2==0) = nan;

ftsize = 15;
%%
g = figure('Units','points','Position',[50,50,750,350]);

subplot(1,2,1)
p1 = pcolor(x,y,20*log10(abs(field_ref_ps-field_ps)));
set(gca, 'Units','normalized','Position',[ 0.075 0.075 0.4 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
plot(x_ps.*w_ps_2,y_ps.*w_ps_2,'--w','LineWidth',1);
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
p2 = pcolor(x,y,20*log10(abs(field_ref_ps_2-field_ps_2)));
set(gca, 'Units','normalized','Position',[ 0.575,0.075 0.4 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
plot(x_ps2.*w_ps_2,y_ps2.*w_ps_2,'--w','LineWidth',1);
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
ylabel(c,'[dB]');
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

set(gcf,'PaperPositionMode','auto');
print -dpng fixed_referencing_circular_ps -r300