clear;
close all
%
addpath(genpath('../Files'));

[pot_bb, pot_path] = read_epspath('general_contour.eps');
siz = norm(diff(pot_bb));
rad_int = meshpath(pot_path, siz/150);
rad_int = translate_mesh(scale_mesh(rad_int, 1/siz*5.25), [1.5 .5]);
rad_ext = rad_int;
rad_ext.Elements = flip_elements(rad_ext.Elements);

Lx = 6.25;
Ly = 4;
dx = 1e-2;
x = (0:dx:Lx);  
y = (0:dx:Ly);
[X,Y] = meshgrid(x,y);

c = 343.1;
f0 = 1e3;
k = 2*pi *f0/c;
xs = [0.4 2.5];
xr = [3.5, 1.5];

R = sqrt((X-xs(1)).^2 + (Y-xs(2)).^2);
field = 1/(4*pi)*exp(-1i*k*R)./R;

%// extract the mesh description matrices    
x2 = rad_int.Nodes(:,2);
y2 = rad_int.Nodes(:,3);
z2 = rad_int.Nodes(:,4);

in = inpolygon(X(:),Y(:),x2,y2);
mask = double(  reshape ( in, length(y), length(x) ) ) ;


fig = figure('Units','points','Position',[200,200,450,250]);
p = axes('Units','normalized','Position',[0.05 0.025 .95 .95 ]);

pcolor(x,y,real(field));
axis equal tight
caxis([-1 1] * 1e-1);

hold on
mask(mask == 0)= nan;
shading interp

f1 = fill([x2', x2(end) Lx Lx 0 0 x2(end) x2(end) ],[y2' Ly Ly 0 0 Ly Ly y2(end)],[1 1 1]*1,'EdgeColor','none','FaceAlpha',0.75);

line(x2, y2, z2, 'Color', 'black','LineStyle','-','LineWidth',2);

pcolor_ax =(gca);
cRange = caxis;
[C,~]= contour( x, y , real(field), '-k');


line( [xs(1),xr(1)],[xs(2),xr(2)],'Color','black', 'LineStyle','--' )


N = 100;
L = [ linspace(xs(1),xr(1),100)' linspace(xs(2),xr(2),100)'  ];
val = zeros(N,1);
ind = zeros(N,1);
for n = 1 : N
    [val(n),ind(n)] = min(sqrt( (x2 - L(n,1)).^2 + (y2 - L(n,2)).^2 )); 
end
[~,i] = min(val);

x0 = [x2(ind(i)),y2(ind(i))];

R0 = sqrt((X-xr(1)).^2 + (Y-xr(2)).^2);
field2 =exp(-1i*k*R0)./R0;
[C2,~]= contour( x, y , real(field2).*mask, '--k');
set(gca,'XTickLabel','','YTickLabel','');


plot(xs(1),xs(2),'.k','MarkerSize',20);
plot(xr(1),xr(2),'.k','MarkerSize',20);
plot(x0(1),x0(2),'.k','MarkerSize',20);


kxn_p = -(x0(1)-xs(1))/sqrt( sum( (x0-xs).^2 ) );
kyn_p = -(x0(2)-xs(2))/sqrt( sum( (x0-xs).^2 ) );
kxn_g = -(x0(1)-xr(1))/sqrt( sum( (x0-xr).^2 ) );
kyn_g = -(x0(2)-xr(2))/sqrt( sum( (x0-xr).^2 ) );
ftsize = 12;


headWidth = 5*1.5;
headLength = 5*1.5;

LineLength = 0.65;
ah3 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth,'LineWidth',1.5);
set(ah3,'parent',gca);
set(ah3,'position',[x0(1), x0(2), LineLength*kxn_p, LineLength*kyn_p]);    

ah4 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth,'LineWidth',1.5);
set(ah4,'parent',gca);
set(ah4,'position',[x0(1), x0(2), LineLength*kxn_g, LineLength*kyn_g]);    

ah5 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth,'LineWidth',1.5);
set(ah5,'parent',gca);
set(ah5,'position',[xr(1), xr(2), -LineLength*kxn_p, -LineLength*kyn_p]);    


set(gca,'box','off')
color = get(fig,'Color');

set(gca,'XColor',color,'YColor',color,'TickDir','out')
set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/High_freq_approximations','KHIE_stat_point' ) ,'-dpng')