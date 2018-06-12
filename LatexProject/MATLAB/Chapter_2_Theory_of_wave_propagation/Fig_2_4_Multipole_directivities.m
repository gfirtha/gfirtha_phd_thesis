clear all
close all

fi = linspace(0,2*pi,50);
theta = linspace(-pi,pi,50);
[Fi,Theta] = meshgrid(fi,theta);

monopole = 1;
dipole = sin(Fi).*cos(Theta);
quadrupole = sin(2*Fi).*cos(Theta);


[X1,Y1,Z1] = sph2cart(Fi, Theta ,abs(monopole));
[X2,Y2,Z2] = sph2cart(Fi, Theta ,abs(dipole));
[X3,Y3,Z3] = sph2cart(Fi, Theta ,abs(quadrupole));
%%
ftsize = 11;
f = figure('Units','points','Position',[200,200,500,180]);
set(f,'defaulttextinterpreter','latex')

pos = [ 0.0    0.15 0.33 .8 
        0.35  0.15 0.3 .85
        0.695   0.15 0.3 .85 ];
        
    
p1 = axes('Units','normalized','Position',pos(1,:));
set(gca,'TickLabelInterpreter', 'tex');
surf(X1,Y1,Z1);
axis equal tight
xlabel('$x$','FontSize',ftsize);
ylabel('$y$','FontSize',ftsize);
zlabel('$z$','FontSize',ftsize);
set(gca,'FontName','Times New Roman');
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');

p2 = axes('Units','normalized','Position',pos(2,:));
set(gca,'TickLabelInterpreter', 'tex');
surf(X2,Y2,Z2);
axis equal tight
xlabel('$x$','FontSize',ftsize);
ylabel('$y$','FontSize',ftsize);
zlabel('$z$','FontSize',ftsize);
set(gca,'FontName','Times New Roman');
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');

p3 = axes('Units','normalized','Position',pos(3,:));
set(gca,'TickLabelInterpreter', 'tex');
surf(X3,Y3,Z3);
axis equal tight
xlabel('$x$','FontSize',ftsize);
ylabel('$y$','FontSize',ftsize);
zlabel('$z$','FontSize',ftsize);
set(gca,'FontName','Times New Roman');
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');


set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');

set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/Basic_acoustics','monopole_dipole' ) ,'-dpng')