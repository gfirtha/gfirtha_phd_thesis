clear
close all

dx = 1e-2;
w = 2*pi*1e3;
c = 343.1;
k = w / c;

xr = [1.25 , 1.35];
xs = [0,-1];
x0 = [(xr(1)-xs(1))/(1+abs(xr(2)/xs(2)))+xs(1), 0];
x = (-2:dx:2)';
y = (-1.5:dx:2)';
[X,Y] = meshgrid(x,y);
r = sqrt( (X-xs(1)).^2 + (Y-xs(2)).^2 );

angle_0 = 1i*(k*sqrt( sum( (x0-xs).^2 ) )+pi/2);
field = exp( -1i*k*r )./r;

angle_1 = -1i*k*sqrt( sum((x0-xs).^2 ) ) ;
angle_2 = -1i*k*sqrt( sum((x0-xr).^2) );
r2 = sqrt( (X-x0(1)).^2+(Y-x0(2)).^2);
field2 = exp(-1i*k*r2 +  angle_1)./r2;
kx0 = (x0-xs)/sqrt( sum( (x0-xs).^2 ) );
ky0 = (x0-xs)/sqrt( sum( (x0-xs).^2 ) );
%%
x_ax = -1.5;
fig = figure('Units','points','Position',[200,200,400,330]);
p = axes('Units','normalized','Position',[0.05 0.025 .95 .95 ]);
p1 = pcolor(x,y,real(field));
set(p1, 'EdgeColor', 'white');
shading interp
axis equal tight
caxis([-1,1]*1e0*1.5);
hold on
fill([x(1), x(end),x(end),x(1)],[y(1),y(1),0,0],[1 1 1]*1,'EdgeColor','none','FaceAlpha',0.75);
fill([x(1), x(end),x(end),x(1)],[0,0,y(end),y(end)],[1 1 1]*1,'EdgeColor','none','FaceAlpha',0);

plot(xs(1),xs(2),'.k','MarkerSize',20)
plot(xr(1),xr(2),'.k','MarkerSize',20)
plot(x0(1),x0(2),'.k','MarkerSize',20)
plot(xr(1),0,'.k','MarkerSize',20)

pcolor_ax =(gca);
cRange = caxis;
[C,~]= contour( x, y , real(field), [ 0 0 ], '-k');
field2(r2>sqrt(sum(( xr - x0 ).^2))+10*dx)=nan;
field2(Y<0) = nan;

[C2,~]= contour( x, y , real(field2), '--k');
set(gca,...
'XTickLabel','','YTickLabel','');
line([xs(1),xr(1)],[xs(2),xr(2)],'Color','black','LineStyle','--')
line([x(1),x(end)],[0,0],'Color','black','LineStyle','-','LineWidth',2)
line([x(1),x(end)],[xr(2),xr(2)],'Color','black','LineStyle','--','LineWidth',1)

line([xr(1),xr(1)],[0,xr(2)],'Color','black','LineStyle',':','LineWidth',1)


headWidth = 5*1.5;
headLength = 5*1.5;
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[x(1)-10*dx, 0, x(end)-x(1)+25*dx, 0]);    

ah2 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah2,'parent',gca);
set(ah2,'position',[x_ax, y(1), 0, y(end)-y(1)+15*dx]); 

kxn_p = -(x0(1)-xs(1))/sqrt( sum( (x0-xs).^2 ) );
kyn_p = -(x0(2)-xs(2))/sqrt( sum( (x0-xs).^2 ) );
kxn_g = -(x0(1)-xr(1))/sqrt( sum( (x0-xr).^2 ) );
kyn_g = -(x0(2)-xr(2))/sqrt( sum( (x0-xr).^2 ) );

ftsize = 12;
LineLength = 0.65;
% ah3 = annotation('arrow',...
%     'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth,'LineWidth',1.5);
% set(ah3,'parent',gca);
% set(ah3,'position',[x0(1), x0(2), LineLength*kxn_p, LineLength*kyn_p]);    
ah4 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth,'LineWidth',1.5);
set(ah4,'parent',gca);
set(ah4,'position',[xr(1), xr(2), LineLength*kxn_g, LineLength*kyn_g]);    
%text(x(end)+0.15,-.1,'$x$',  'Interpreter', 'LaTex' , 'FontSize', ftsize)
%text(x_ax-0.17,y(end)+0.15,'$y$',  'Interpreter', 'LaTex' , 'FontSize', ftsize)

set(gca,'box','off')
color = get(fig,'Color');
xlim([x(1),x(end)+20*dx])
ylim([y(1),y(end)+20*dx])
set(gca,'XColor',color,'YColor',color,'TickDir','out')
set(gca,'FontName','Times New Roman','FontSize',12);
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/SFS_theory','explicit_sol_stationary_point' ) ,'-dpng')