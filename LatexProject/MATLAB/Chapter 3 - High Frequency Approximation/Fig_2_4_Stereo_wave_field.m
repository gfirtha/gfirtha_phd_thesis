clear
close all
addpath(genpath('../Files/arrow3'));
dx = 1e-2;
x = (-2:dx:2)';
y = (-1:dx:2.5)';

w = 0.8e3*2*pi;
c = 343.1;
k = w / c;

R0 = 2.5;
fi0 = 30;

x1 = R0*[-sind(fi0) cosd(fi0)];
x2 = R0*[ sind(fi0) cosd(fi0)];

fi = 10;
C = ( tan(fi*pi/180)/tand(fi0) );
A1 = 1;
A2 = A1*(1-C)/(C+1);
%A1 = 1;A2 = ;
%
Avec = [A1 A2];
A1 = A1/norm(Avec);  A2 = A2/norm(Avec);

[X,Y] = meshgrid(x,y);
r1 = sqrt((X-x1(1)).^2 + (Y-x1(2)).^2);
r2 = sqrt((X-x2(1)).^2 + (Y-x2(2)).^2);

r1px = (X-x1(1))./r1;                       r2px = (X-x2(1))./r2;
kx = (k*(r1px.*A1^2.*r2.^2+ r2px.*A2^2.*r1.^2 + A1*A2*r1.*r2.*(r1px+r2px).*cos(k*(r1-r2)))...-
    - A1*A2*(r1px.*r2-r2px.*r1).*sin(k*(r1-r2)))./((A1*r2).^2+(A2*r1).^2+2*A1*A2.*r1.*r2.*cos(k*(r1-r2)));
r1py = (Y-x1(2))./r1;                       r2py = (Y-x2(2))./r2;
ky = (k*(r1py.*A1^2.*r2.^2+ r2py.*A2^2.*r1.^2 + A1*A2*r1.*r2.*(r1py+r2py).*cos(k*(r1-r2)))...-
    - A1*A2*(r1py.*r2-r2py.*r1).*sin(k*(r1-r2)))./((A1*r2).^2+(A2*r1).^2+2*A1*A2.*r1.*r2.*cos(k*(r1-r2)));
K = sqrt(kx.^2+ky.^2); 
dir = atan(kx./ky);
kxn = kx./K;                                 kyn = ky./K;

[kxx, ~] = gradient(kx,dx);
[~, kyy] = gradient(ky,dx);
curv = kxx + kyy;

field = 1/(4*pi)*(A1*exp(-1i*k*r1)./r1 + A2*exp(-1i*k*r2)./r2);

%%
ftsize= 12;
f = figure('Units','points','Position',[200,200,520,300]);
pos = [ 0.02 0.15 0.6 .75
        0.71 0.73 0.28 .16
        0.71 0.45 0.28 .16
        0.71 0.14 0.28 .16];

p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field));
shading interp
axis equal tight
caxis(1*[-1,1]*1e-1);
hold on
pcolor_ax =(gca);
cRange = caxis; 
[C,~] = contour( x, y , real(field), '-k');
hLines = findobj(gca, 'type', 'line'); % find all the separate lines on contour plot.
set(hLines, 'LineWidth', 1); % and set their width.
caxis(cRange); 
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
plot(y(1:5:end)*0,y(1:5:end),'--k');
% S = contourdata(C);
% %N=6,10
% N = 13;
% x0 = S(N).xdata;
% y0 = S(N).ydata;
% plot(x0,y0,'--k');
% 
x0 = 0*y;
y0 = y;
kx0 = kx(:,x == 0);
ky0 = ky(:,x == 0);
kxn0 = kxn(:,x == 0);
kyn0 = kyn(:,x == 0);
curv0 = curv(:,x == 0);

headWidth = 5;
headLength = 5;
LineLength = 0.5;
for i = 1:15:length(kx0)
    %if K(i)< k*1.5
        ah = annotation('arrow',...
            'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[x0(i) y0(i) LineLength*kxn0(i) LineLength*kyn0(i)]);
  %  end
end
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
xtp1 = p1.XTickLabel;

p2 = axes('Units','normalized','Position',pos(2,:));
plot(y,kx0/k)
xlabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$k_x^P(x=0,y)/k$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2 );
xlim([y(1),y(end)])
grid on
set(gca,'FontName','Times New Roman');
set (gca,'Xdir','reverse')
hold on
yl = ylim;
plot(x1(2)+0*linspace(yl(1),yl(2),20),linspace(yl(1),yl(2),20),'.k','MarkerSize',0.25)

p3 = axes('Units','normalized','Position',pos(3,:));
plot(y,-ky0/k)
xlabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$k_y^P(x=0,y)/k$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2 );
xlim([y(1),y(end)])
grid on
set(gca,'FontName','Times New Roman');
set (gca,'Xdir','reverse')
hold on
yl = ylim;
plot(x1(2)+0*linspace(yl(1),yl(2),20),linspace(yl(1),yl(2),20),'.k','MarkerSize',0.25)

p4 = axes('Units','normalized','Position',pos(4,:));
plot(y,sqrt(kx0.^2+ky0.^2)/k)
xlabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$\left| \mathbf{k}^P(x=0,y)/k \right|$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize-2 );
xlim([y(1),y(end)])
grid on
set(gca,'FontName','Times New Roman');
set (gca,'Xdir','reverse')
hold on
yl = ylim;
plot(x1(2)+0*linspace(yl(1),yl(2),20),linspace(yl(1),yl(2),20),'.k','MarkerSize',0.25)

allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);
set(p1,'XTickLabel',xtp1)
set(gcf,'PaperPositionMode','auto');
print( fullfile( '../..','Figures/High_freq_approximations','stereophony' ) ,'-dpng')