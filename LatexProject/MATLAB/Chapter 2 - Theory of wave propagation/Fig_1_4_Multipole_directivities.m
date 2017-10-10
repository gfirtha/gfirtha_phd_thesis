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
ftsize = 14;
f = figure('Units','points','Position',[200,200,750,240]);

subplot(1,3,1)
set(gca,'TickLabelInterpreter', 'tex');
set(gca, 'Units','normalized','Position',[ 0.012 0.15 0.33 .8 ]);
surf(X1,Y1,Z1);
axis equal tight
xlabel('$x$','Interpreter','latex','FontSize',ftsize);
ylabel('$y$','Interpreter','latex','FontSize',ftsize);
zlabel('$z$','Interpreter','latex','FontSize',ftsize);
set(gca,'FontName','Times New Roman');


subplot(1,3,2)
set(gca,'TickLabelInterpreter', 'tex');
set(gca, 'Units','normalized','Position',[ 0.33 0.15 0.33 .8 ]);
surf(X2,Y2,Z2);
axis equal tight
xlabel('$x$','Interpreter','latex','FontSize',ftsize);
ylabel('$y$','Interpreter','latex','FontSize',ftsize);
zlabel('$z$','Interpreter','latex','FontSize',ftsize);
set(gca,'FontName','Times New Roman');


subplot(1,3,3)
set(gca,'TickLabelInterpreter', 'tex');
set(gca, 'Units','normalized','Position',[ 0.685 0.15 0.33 .8 ]);
surf(X3,Y3,Z3);
axis equal tight
xlabel('$x$','Interpreter','latex','FontSize',ftsize);
ylabel('$y$','Interpreter','latex','FontSize',ftsize);
zlabel('$z$','Interpreter','latex','FontSize',ftsize);
set(gca,'FontName','Times New Roman');


set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
%set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);


set(gcf,'PaperPositionMode','auto');
print( fullfile( '../..','Figures/Basic_acoustics','monopole_dipole' ) ,'-dpng')