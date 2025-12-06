clear all
close all

dx = 0.025;
x0 = (-70:dx:70)';
xs = [0,-1];
alpha = 60*pi/180;
y = (dx:dx:3)';
yref = 2;

[X,Y] = meshgrid(x0,y);
f = 1.5e3;
speed_of_sound = 343.1;
omega = 2*pi*f;
k = omega/speed_of_sound;

r_s = sqrt( X.^2 + ( Y-xs(2) ).^2 );   % virtual source to reference line
r_0 = sqrt( X.^2 + Y.^2 );              
r_2 = sqrt( x0.^2 + xs(2)^2 );             % virtual source to ssd distance

P_ref_ls =        - 1i/4 * besselh(0,2,k*r_s);
P_ref_pw =        exp( -1i*k*(cos(alpha)*X + sin(alpha)*Y) );
G_ssd_to_yref  =   1/(4*pi) * exp(-1i*k*r_0)./r_0;

ys = abs(xs(2));
r0 = sqrt(x0.^2 + ys.^2);

% Reference on a line
d_ls = yref*r0./ys;
d_pw = yref/sin(alpha);

% Reference on a circle
%Rref = ys + yref;
%d_ls = Rref-r0; 
%d_pw = -cos(alpha)*x0 + sqrt( (Rref-ys)^2 - (sin(alpha)*x0).^2 );

% Fixed reference distance (Ahrens)
%d_ls = yref;
%d_pw = yref;

% Reference point (Spors: WFS revisited)
%xref = [0 yref];
%d_ls = sqrt( x0.^2 + yref^2 );
%d_pw = sqrt( x0.^2 + yref^2 );

D_ls = sqrt(2*pi/(1i*k))*sqrt(d_ls)*0.5*1i*k*xs(2)./r_2.*besselh(1,2,k*r_2);
D_pw = sqrt(8*pi/(1i*k))*sqrt(d_pw)*1i*k*sin(alpha).*exp(-1i*k*cos(alpha)*x0);

P_synth_ls = zeros(size(X));
P_synth_pw = zeros(size(X));

wb = waitbar(0,'Calculating radiated field');
for n = 1:length(y)
    waitbar(n/length(y),wb);
    P_synth_ls(n,:)= conv(D_ls,G_ssd_to_yref(n,:),'same')*dx;
    P_synth_pw(n,:)= conv(D_pw,G_ssd_to_yref(n,:),'same')*dx;
end
close(wb)

x_c_ls = x0 + d_ls.*x0./r0;
y_c_ls = d_ls.*ys./r0;

x_c_pw = x0 + d_pw.*cos(alpha);
y_c_pw = d_pw.*sin(alpha);

%%
f = figure('Units','points','Position',[200,200,730,200]);
set(gcf,'Units','normalized');
ftsize = 16;

x_lim = 3;
i = find(abs(x0)<= x_lim);
subplot(1,2,1)
p1 = pcolor(x0(i),y,real(P_synth_pw(:,i)));shading interp;axis equal tight
set(gca, 'Units','normalized','Position',[ 0.06 0.075 .4 .9 ]);
caxis([-2,2])
ylim([y(1),y(end)])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$\mathcal{R}(P_{\mathrm{synth,pw}}(\mathbf{x},\omega))$', 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-6);

subplot(1,2,2)
p2 = pcolor(x0(i),y,real(P_synth_ls(:,i)));shading interp;axis equal tight
set(gca, 'Units','normalized','Position',[ 0.56 0.075 .4 .9 ]);
caxis([-.05,.05])

ylim([y(1),y(end)])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$\mathcal{R}(P_{\mathrm{synth,ls}}(\mathbf{x},\omega))$', 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-6);
set(gcf,'PaperPositionMode','auto');
%print -dpng real_part -r300
%%
f = figure('Units','points','Position',[200,50,600,370]);

x_lim = 6;
i = find(abs(x0)<= x_lim);
set(gcf,'Units','normalized');
ftsize = 14;
subplot(2,1,1)
p1 = pcolor(x0(i),y,20*log10(abs(P_ref_pw(:,i)-P_synth_pw(:,i))));shading interp;axis equal tight
set(gca, 'Units','normalized','Position',[ 0.065 0.575 .9 .4 ]);
c = colorbar;
ylabel(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize)
caxis([-30,20])
xlim([-6,6])
hold on
plot(x_c_pw,0*x_c_pw+y_c_pw,'--w','LineWidth',1)
ylim([y(1),y(end)])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,pw}}(\mathbf{x},\omega)-P_{\mathrm{ref,pw}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);

subplot(2,1,2)
p2 = pcolor(x0(i),y,20*log10(abs(P_ref_ls(:,i)-P_synth_ls(:,i))));shading interp;axis equal tight
set(gca, 'Units','normalized','Position',[ 0.065 0.075 .9 .4 ]);
c = colorbar;
ylabel(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize)
caxis([-70,0])
hold on
plot(x_c_ls,y_c_ls,'--w','LineWidth',1)
xlim([-6,6])
ylim([y(1),y(end)])
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{synth,ls}}(\mathbf{x},\omega)-P_{\mathrm{ref,ls}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);
set(gcf,'PaperPositionMode','auto');
%print -dpng circle_referencing -r300