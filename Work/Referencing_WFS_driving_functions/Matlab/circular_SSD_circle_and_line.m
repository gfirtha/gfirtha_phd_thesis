clear all
close all

R0 = 1.5;
fi = linspace(0,2*pi,1000)';
dfi = fi(2)-fi(1);
x_ssd = R0*cos(fi);
y_ssd = R0*sin(fi);
n_ssd = [ -cos(fi) -sin(fi) ];

x = linspace(-2,2,200);
y = linspace(-2,2,200);
[X,Y] = meshgrid(x,y);

k = 2e3*2*pi/343.1;

% Reference a plane wave to a circle
alpha_pw = 0*pi/180;
a0 = double((n_ssd * [cos(alpha_pw); sin(alpha_pw)])>0);
r_ref = 1;
d_pw =  R0*cos(abs(fi-alpha_pw-pi)) -...
    (sqrt( r_ref.^2 - (sin(pi-fi+alpha_pw)*R0).^2 ) );
alpha_fi = cos(fi)*cos(alpha_pw) + sin(fi)*sin(alpha_pw);
D_pw = -(alpha_fi).*sqrt(8*pi*1i*k).*sqrt(d_pw).*...
    exp(-1i*k*R0*alpha_fi).*a0;

% Reference a point source to a line
x_s = -3;
y_s =  0;
v = [(x_ssd - x_s) (y_ssd - y_s)];
a_ps = dot(v,n_ssd,2)>0; 
r_ps = sqrt( dot(v,v,2) );
v_n = bsxfun( @times, v, 1./r_ps );
alpha = acos(v_n(:,1)).*sign(v_n(:,2));

x0 = -0.5;
d_ps = -r_ps.*( (x0+R0*cos(fi-pi))./(x_s + R0*cos(fi-pi)) );

D_ps = -sqrt(1i*k/(2*pi)).*(cos(fi).*v(:,1)+sin(fi).*v(:,2)).*...
        sqrt(d_ps.*r_ps./(d_ps+r_ps)).*a_ps.*exp(-1i*k*r_ps)./r_ps.^2;
%%
field_pw = zeros(size(X));
field_ps = zeros(size(X));

for n = 1:length(x_ssd)
   
    r_ssd = sqrt(( X - x_ssd(n) ).^2 + ( Y - y_ssd(n) ).^2 );
    
    field_pw = field_pw + R0/(4*pi)*D_pw(n)*exp(-1i*k*r_ssd)./r_ssd*dfi;
                                        
    field_ps = field_ps + R0/(4*pi)*D_ps(n)*exp(-1i*k*r_ssd)./r_ssd*dfi;
    
end

field_pw_ref = exp(-1i*k*(cos(alpha_pw)*X+sin(alpha_pw)*Y));

R = sqrt( (X-x_s).^2 + (Y-y_s).^2 );
field_ps_ref = 1/(4*pi)*exp(-1i*k*R)./R;
%%
x_pw = x_ssd + d_pw*cos(alpha_pw);
y_pw = y_ssd + d_pw*sin(alpha_pw);

x_ps = x_s + (d_ps+r_ps).*cos(alpha);
y_ps = y_s + (d_ps+r_ps).*sin(alpha);
%%
ftsize = 15; 
f = figure('Units','points','Position',[50,50,750,350]);
subplot(1,2,1)
p1 = pcolor(x,y,20*log10(abs(field_pw_ref-field_pw)));
set(gca, 'Units','normalized','Position',[ 0.075 0.075 0.4 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
plot(x_pw,y_pw,'--w','LineWidth',1)
caxis([-30,20])
c = colorbar;
ylabel(c,'[dB]')
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,pw}}(\mathbf{x},\omega)-P_{\mathrm{ref,pw}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
c = colorbar;
ylabel(c,'[dB]')
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

subplot(1,2,2)
p2 = pcolor(x,y,20*log10(abs(field_ps_ref-field_ps)));
set(gca, 'Units','normalized','Position',[ 0.575,0.075 0.4 .9 ]);
shading interp
axis equal tight
hold on
plot(x_ssd,y_ssd,'k','LineWidth',1)
plot(x_ps,y_ps,'--w','LineWidth',1)
caxis([-60,-10])
c = colorbar;
ylabel(c,'[dB]')
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,ps}}(\mathbf{x},\omega)-P_{\mathrm{ref,ps}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
c = colorbar;
ylabel(c,'[dB]')
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

set(gcf,'PaperPositionMode','auto');
print -dpng pw_to_circle_ps_to_line -r300