clear
close all

phi = linspace(0,2*pi,50);
el = linspace(-pi/2,pi/2,50);
[Phi,El] = meshgrid(phi, el);

monopole = ones(size(Phi));
dipole = sin(Phi).*cos(El);
quadrupole = sin(2*Phi).*cos(El);


[X1,Y1,Z1] = sph2cart(Phi, El ,abs(monopole));
[X2,Y2,Z2] = sph2cart(Phi, El ,abs(dipole));
[X3,Y3,Z3] = sph2cart(Phi, El ,abs(quadrupole));
%%
ftsize = 11;
f = figure('Units','points','Position',[200,200,500,180]);
set(f,'defaulttextinterpreter','latex')

pos = [ 0.0     0.15 0.33 .8 
        0.35    0.15 0.3  .85
        0.695   0.15 0.3  .85 ];
        
    
p1 = axes('Units','normalized','Position',pos(1,:));
set(gca,'TickLabelInterpreter', 'tex');
surf(X1,Y1,Z1,monopole);
axis equal tight
xlabel('$x$','FontSize',ftsize);
ylabel('$y$','FontSize',ftsize);
zlabel('$z$','FontSize',ftsize);
set(gca,'FontName','Times New Roman');
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');
caxis([-1,1])

p2 = axes('Units','normalized','Position',pos(2,:));
set(gca,'TickLabelInterpreter', 'tex');
surf(X2,Y2,Z2,dipole);
axis equal tight
xlabel('$x$','FontSize',ftsize);
ylabel('$y$','FontSize',ftsize);
zlabel('$z$','FontSize',ftsize);
set(gca,'FontName','Times New Roman');
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');
caxis([-1,1])

p3 = axes('Units','normalized','Position',pos(3,:));
set(gca,'TickLabelInterpreter', 'tex');
surf(X3,Y3,Z3,quadrupole);
axis equal tight
xlabel('$x$','FontSize',ftsize);
ylabel('$y$','FontSize',ftsize);
zlabel('$z$','FontSize',ftsize);
set(gca,'FontName','Times New Roman');
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');
caxis([-1,1])

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/Basic_acoustics','monopole_dipole' ) ,'-dpng')