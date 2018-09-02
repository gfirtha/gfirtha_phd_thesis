clear
close all

c = 343.1;      % Speed of sound [m/s]
f = 1000;       % Plane wave frequency
w = 2*pi*f;     % Angular frequency
k = w/c;        % Acoustic wave number
lambda = c/f;

kx1 = k*1/sqrt(2);
kz1 = 0;
kx2 = k*1.01;
kz2 = 0;

ky1 = sqrt( k^2 - kx1^2 - kz1^2 );
ky2 = sqrt( k^2 - kx2^2 - kz2^2 );

dx = 0.01;
x = (-3.1*lambda:dx:3.1*lambda)';
y = x+3.1*lambda;
[X,Y] = meshgrid(x,y);

field_prop = exp( 1i*( kx1*X + ky1*Y  ) );
field_evan = exp( 1i*( kx2*X + ky2*Y  ) );
%%
ftsize = 11/0.9;
f = figure('Units','points','Position',[200,200,500,230]);
set(f,'defaulttextinterpreter','latex')
set(gcf,'Units','normalized');

subplot(1,2,1)
p1 = pcolor(x/lambda,y/lambda,real(field_prop));
set(gca, 'Units','normalized','Position',[ 0.09 0.11 0.37 .9 ]);
caxis([-1.5,1.5])
shading interp
xlabel( '$x/\lambda \rightarrow$ []' , 'FontSize', ftsize );
ylabel( '$y/\lambda \rightarrow$ []' , 'FontSize', ftsize );
axis equal tight
hold on
contour( x/lambda, y/lambda, real(field_prop),[0 0], '-k');

c = [-.25 .75];
xo = [0,kx1]/k;
yo = [0,ky1]/k;

headWidth = 5;
headLength = 5;
LineLength = 0.5;
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[c(1) c(2) kx1/k  ky1/k ]/lambda);
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[c(1) c(2) kx1/k 0]/lambda);
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[c(1) c(2) 0 ky1/k]/lambda);
plot(c(1)/lambda,c(2)/lambda, 'ko', 'MarkerSize',2,'LineWidth',2);
text( 0.28/lambda, 1.45/lambda , 0.5 ,'$\mathbf{k}$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
text( 0.35/lambda, 0.85/lambda , 0.5 ,'$k_x$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
text( -0.2/lambda, 1.4/lambda , 0.52 ,'$k_y$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
text(  -0.1/lambda, 0.85/lambda , 0.52 ,'$\varphi$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );

t = (-0.75:0.01:0.75)';
plot(c(1)/lambda-t,c(2)/lambda+t,'--k');
fi = (0:0.1:45)*pi/180;
plot( 0.3*cos(fi)+c(1)/lambda, 0.3*sin(fi)+c(2)/lambda , 'k' );
set(gca,'FontName','Times New Roman');


subplot(1,2,2)
p2 = pcolor(x/lambda,y/lambda,real(field_evan));
set(gca, 'Units','normalized','Position',[ 0.61 0.11 0.37 .9 ]);
caxis(.5*[-1,1])
shading interp
xlabel( '$x/\lambda \rightarrow []$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$y/\lambda \rightarrow []$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
axis equal tight
hold on
pcolor_ax =(gca);
cRange = caxis;
[C,~]= contour( x/lambda, y/lambda, real(field_evan),[0 0], '-k');

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/Basic_acoustics','plane_wave_illustration' ) ,'-dpng')