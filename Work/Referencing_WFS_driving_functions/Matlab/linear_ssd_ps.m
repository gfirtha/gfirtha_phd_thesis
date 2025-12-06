clear all
close all

dx = 0.025;
x0 = (-15:dx:15)';
xs = [0,-3];

y = (dx:dx:4)';
yref =  1.5;

[X,Y] = meshgrid(x0,y);
f = 1.5e3;
speed_of_sound = 343.1;
omega = 2*pi*f;
k = omega/speed_of_sound;

r_s = sqrt( X.^2 + ( Y-xs(2) ).^2 );   % virtual source to reference line
dr  = sqrt( X.^2 + Y.^2 );              

Gref =            1/(4*pi)*exp(-1i*k*r_s)./r_s;
G_ssd_to_yref  =  1/(4*pi)*exp(-1i*k* dr)./dr;

ys = abs(xs(2));
r0 = sqrt( x0.^2 + ys.^2);

% Fixed reference (Ahrens)
d_fixed = yref;

%  Reference on a line
d_line = r0*yref/(ys+yref);

% Reference on a circle
Rref = yref+ys;
d_circle = r0.*(Rref-r0)/Rref;

% Spors Revisited referencing
%xref = [0 yref];
%d = sqrt( x0.^2 + yref^2 );

dP = -sqrt(1i*k/(2*pi)).*xs(2).*exp(-1i*k*r0)./r0.^2;
D_fixed =  sqrt(d_fixed).*dP;
D_line  =  sqrt(d_line).*dP;
D_circle = sqrt(d_circle).*dP;

P_fixed = zeros(size(X));
P_line = zeros(size(X));
P_circle = zeros(size(X));
wb = waitbar(0,'Calculating ');
for n = 1:length(y)
    waitbar(n/length(y),wb);
    P_fixed(n,:) = conv(D_fixed ,G_ssd_to_yref(n,:),'same')*dx;
    P_line(n,:)  = conv(D_line  ,G_ssd_to_yref(n,:),'same')*dx;
    P_circle(n,:)= conv(D_circle,G_ssd_to_yref(n,:),'same')*dx;
end
close(wb)
%%
d_ps_1 = d_fixed.*r0./(r0-d_fixed);
x_c_fixed = x0 + d_ps_1.*x0./r0;
y_c_fixed = d_ps_1.*ys./r0;

d_ps_2 = d_line.*r0./(r0-d_line);
x_c_line = x0 + d_ps_2.*x0./r0;
y_c_line = d_ps_2.*ys./r0;

d_ps_3 = d_circle.*r0./(r0-d_circle);
x_c_circle = x0 + d_ps_3.*x0./r0;
y_c_circle = d_ps_3.*ys./r0;
%%

f = figure('Units','points','Position',[100,30,500,500]);
set(gcf,'Units','normalized');
ftsize = 12;
xind = find(abs(x0)<=6) ;

subplot(3,1,1)
sp1 = pcolor( x0(xind), y, 20*log10( abs( Gref(:,xind) - P_fixed(:,xind) ) ) );
set(gca, 'Units','normalized','Position',[ 0.1 0.69 .85 .25 ]);
hold on
shading interp;axis equal tight
c = colorbar;
ylabel(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize)
caxis([-80,10])
hold on
plot(x_c_fixed,y_c_fixed,'--w','LineWidth',1)
xlim([-6,6]);
ylim([y(1),y(end)]);
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{fixed,ps}}(\mathbf{x},\omega)-P_{\mathrm{ref,ps}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);


subplot(3,1,2)
sp2 = pcolor( x0(xind), y, 20*log10( abs( Gref(:,xind) - P_line(:,xind) ) ) );
set(gca, 'Units','normalized','Position',[ 0.1 0.37 .85 .22 ]);
hold on
shading interp;axis equal tight
c = colorbar;
ylabel(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize)
caxis([-70,10])
hold on
plot(x_c_line,y_c_line,'--w','LineWidth',1)
xlim([-6,6]);
ylim([y(1),3]);
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{linear,ps}}(\mathbf{x},\omega)-P_{\mathrm{ref,ps}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);


subplot(3,1,3)
sp3 = pcolor( x0(xind), y, 20*log10( abs( Gref(:,xind) - P_circle(:,xind) ) ) );
set(gca, 'Units','normalized','Position',[ 0.1 0.07 .85 .22 ]);
hold on
shading interp;axis equal tight
c = colorbar;
ylabel(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize)
caxis([-70,10])
hold on
plot(x_c_circle,y_c_circle,'--w','LineWidth',1)
xlim([-6,6]);
ylim([y(1),3]);
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
title('$20\mathrm{log}10\left( P_{\mathrm{circle,ps}}(\mathbf{x},\omega)-P_{\mathrm{ref,ps}}(\mathbf{x},\omega) \right)$'...
                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);

set(gcf,'PaperPositionMode','auto');
%print -dpng point_source_referencing -r300