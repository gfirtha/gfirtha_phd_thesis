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
alpha = 0*pi/180;
%
w_pw = double((n_ssd * [cos(alpha); sin(alpha)])>0);

% Generalized 2-1/2D WFS
d_fixed = 0.75;

alpha_fi = cos(fi)*cos(alpha) + sin(fi)*sin(alpha);
D_pw = -(alpha_fi).*sqrt(8*pi*1i*k).*sqrt(d_fixed).*...
    exp(-1i*k*R0*alpha_fi).*w_pw;


x_s = -3.0;
y_s =  0;
v = [(x_ssd - x_s) (y_ssd - y_s)]; 
r0 = sqrt(dot(v,v,2));

w_ls = dot(v,n_ssd,2)>0;
D_ls = sqrt(pi*1i*k*d_fixed/2).*(cos(fi).*v(:,1)+sin(fi).*v(:,2)) ./r0.*besselh(1,2,k*r0).*w_ls;

D_ps = -sqrt(1i*k/(2*pi)).*(cos(fi).*v(:,1)+sin(fi).*v(:,2)).*...
        sqrt(d_fixed).*exp(-1i*k*r0)./r0.^2.*w_ls;
%%
field_pw  = zeros(size(X));
field_ls   = zeros(size(X));
field_ps   = zeros(size(X));
for n = 1:length(x_ssd)
    r_ssd = sqrt(( X - x_ssd(n) ).^2 + ( Y - y_ssd(n) ).^2 );
    
    field_pw  = field_pw  + R0/(4*pi)*D_pw(n)...
                                            *exp(-1i*k*r_ssd)./r_ssd*dfi;
    field_ls  = field_ls  + R0/(4*pi)*D_ls(n)...
                                            *exp(-1i*k*r_ssd)./r_ssd*dfi;
    field_ps  = field_ps  + R0/(4*pi)*D_ps(n)...
                                            *exp(-1i*k*r_ssd)./r_ssd*dfi;
    
end
field_ref_pw = exp(-1i*k*(cos(alpha)*X+sin(alpha)*Y));
field_ref_ls = -1i/4*besselh(0,2,k*sqrt( (X-x_s).^2 + (Y-y_s).^2 ));

R = sqrt( (X-x_s).^2 + (Y-y_s).^2 );
field_ref_ps = 1/(4*pi)*exp(-1i*k*R)./R;
%%
x_pw  = x_ssd + d_fixed *cos(alpha);
y_pw  = y_ssd + d_fixed *sin(alpha);
%
v_n = bsxfun( @times, v, 1./r0);
alpha_0 = acos(v_n(:,1)).*sign(v_n(:,2));
x_ls = x_s + (d_fixed+r0).*cos(alpha_0);
y_ls = y_s + (d_fixed+r0).*sin(alpha_0);
%
x_ps = x_s + (r0./(r0-d_fixed)).*cos(alpha_0).*r0;
y_ps = y_s + (r0./(r0-d_fixed)).*sin(alpha_0).*r0;

w_pw(w_pw==0) = nan;
ftsize = 15;
%%

f = figure('Units','points','Position',[50,50,750,350]);

subplot(1,2,1)
p1 = pcolor(x,y,real(field_pw));
set(gca, 'Units','normalized','Position',[ 0.075 0.075 0.375 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
caxis([-1,1])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$\mathcal{R}(P_{\mathrm{synth,pw}}(\mathbf{x},\omega))$', 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-6);

subplot(1,2,2)
p2 = pcolor(x,y,real(field_ls));
set(gca, 'Units','normalized','Position',[ 0.575,0.075 0.375 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
caxis([-.025,.025]);
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$\mathcal{R}(P_{\mathrm{synth,ls}}(\mathbf{x},\omega))$', 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);

set(gcf,'PaperPositionMode','auto');
%print -dpng real_part_circular -r300
%%
g = figure('Units','points','Position',[50,50,750,350]);

subplot(1,2,1)
p1 = pcolor(x,y,20*log10(abs(field_ref_pw-field_pw)));
set(gca, 'Units','normalized','Position',[ 0.075 0.075 0.4 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
plot(x_pw.*w_pw,y_pw.*w_pw,'--w','LineWidth',1);
caxis([-35,20])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,pw}}(\mathbf{x},\omega)-P_{\mathrm{ref,pw}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(g,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);
c = colorbar;
ylabel(c,'[dB]')

subplot(1,2,2)
p1 = pcolor(x,y,20*log10(abs(field_ref_ls-field_ls)));
set(gca, 'Units','normalized','Position',[ 0.575,0.075 0.4 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
plot(x_ls.*w_pw,y_ls.*w_pw,'--w','LineWidth',1);
caxis([-70,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,ls}}(\mathbf{x},\omega)-P_{\mathrm{ref,ls}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(g,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);
c = colorbar;
ylabel(c,'[dB]');
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

set(gcf,'PaperPositionMode','auto');
%print -dpng fixed_referencing_circular -r300
%%

%g = figure('Units','points','Position',[50,50,750,350]);
figure
p1 = pcolor(x,y,20*log10(abs(field_ref_ps-field_ps)));
%set(gca, 'Units','normalized','Position',[ 0.575,0.075 0.4 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
plot(x_ps.*w_pw,y_ps.*w_pw,'--w','LineWidth',1);
caxis([-70,-10])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,ls}}(\mathbf{x},\omega)-P_{\mathrm{ref,ls}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(g,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);
c = colorbar;
ylabel(c,'[dB]');
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

set(gcf,'PaperPositionMode','auto');