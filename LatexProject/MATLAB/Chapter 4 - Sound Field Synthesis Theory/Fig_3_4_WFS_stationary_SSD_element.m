clear;
close all
addpath(genpath('../Files'));
%
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

R = sqrt((X-xs(1)).^2 + (Y-xs(2)).^2);
field = 1/(4*pi)*exp(-1i*k*R)./R;

%// extract the mesh description matrices    
x_ssd = rad_int.Nodes(:,2);
y_ssd = rad_int.Nodes(:,3);
z_ssd = rad_int.Nodes(:,4);

x_ref = (x_ssd-mean(x_ssd))*0.5 + mean(x_ssd);
y_ref = (y_ssd-mean(y_ssd))*0.5 + mean(y_ssd);

Diff_refX = x_ref(1:end-1)-x_ref(2:end); Diff_refY = y_ref(1:end-1)-y_ref(2:end);
D_ref = [ Diff_refX Diff_refY]; D_ref = bsxfun(@times, D_ref , 1./sqrt(Diff_refX.^2+Diff_refY.^2));
Nref_cont = ( [0 1; -1 0] *D_ref')'; 

kx = bsxfun ( @times , [ x_ref-xs(1) y_ref-xs(2) ] , 1./sqrt( (x_ref-xs(1)).^2 +  (y_ref-xs(2)).^2 ) );
mask =   ( sum( kx .*[Nref_cont; 0 0] , 2 ) >= 0 ) ;
cont2 = [x_ref(mask) y_ref(mask)];
cont2 = cont2(1:end-1,:);
%%

xr = [x_ref(70), y_ref(70)];
in = inpolygon(X(:),Y(:),x_ssd,y_ssd);
mask = double(  reshape ( in, length(y), length(x) ) ) ;


fig = figure('Units','points','Position',[200,200,450,250]);
p = axes('Units','normalized','Position',[0.05 0.025 .95 .95 ]);

pcolor(x,y,real(field));
axis equal tight
caxis([-1 1] * 1e-1);

hold on
mask(mask == 0)= nan;
shading interp

f1 = fill([x_ssd', x_ssd(end) Lx Lx 0 0 x_ssd(end) x_ssd(end) ],[y_ssd' Ly Ly 0 0 Ly Ly y_ssd(end)],[1 1 1]*1,'EdgeColor','none','FaceAlpha',0.75);

line(x_ssd, y_ssd, z_ssd, 'Color', 'black','LineStyle','-','LineWidth',2);

line(x_ref, y_ref, z_ssd, 'Color', 'black','LineStyle','--','LineWidth',1.5);
line(cont2(:,1), cont2(:,2), zeros(length(cont2)), 'Color', 'black','LineStyle','-','LineWidth',2);

pcolor_ax =(gca);
cRange = caxis;
[C,~]= contour( x, y , real(field), '-k');


line( [xs(1),xr(1)],[xs(2),xr(2)],'Color','black', 'LineStyle','--','LineWidth',1.5 )


N = 100;
L = [ linspace(xs(1),xr(1),100)' linspace(xs(2),xr(2),100)'  ];

val = zeros(N,1);
ind = zeros(N,1);
for n = 1 : N
    [val(n),ind(n)] = min(sqrt( (x_ssd - L(n,1)).^2 + (y_ssd - L(n,2)).^2 )); 
end
[~,i] = min(val);

x0 = [x_ssd(ind(i)),y_ssd(ind(i))];

R0 = sqrt((X-x0(1)).^2 + (Y-x0(2)).^2);

phase0 = -1i*k* sqrt((x0(1)-xs(1)).^2 + (x0(2)-xs(2)).^2);
field2 =exp(-1i*k*R0+phase0)./R0;
[C2,~]= contour( x, y , real(field2).*mask, '--k');
set(gca,...
'XTickLabel','','YTickLabel','');



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
ah4 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth,'LineWidth',1.5);
set(ah4,'parent',gca);
set(ah4,'position',[x0(1), x0(2), LineLength*kxn_g, LineLength*kyn_g]);    

%%
No_ls = 75;
Dx = x_ssd(2:end)-x_ssd(1:end-1);
Dy = y_ssd(2:end)-y_ssd(1:end-1);
ArcLength = sqrt( Dx.^2 + Dy.^2);
CurveLength = cumsum ( ArcLength );

rep_elem = find(( Dx == 0).*( Dy == 0));
x_ssd_temp = x_ssd;
x_ssd_temp( rep_elem ) = [];
y_ssd_temp = y_ssd;
y_ssd_temp( rep_elem )  = [];
CurveLength( rep_elem )  = [];

x_new_nodes = interp1( [0;CurveLength], x_ssd_temp, linspace(0,CurveLength(end),No_ls+1) )';
y_new_nodes = interp1( [0;CurveLength], y_ssd_temp, linspace(0,CurveLength(end),No_ls+1) )';

x_ls = [(x_new_nodes(1:end-1)+x_new_nodes(2:end))/2 (y_new_nodes(1:end-1)+y_new_nodes(2:end))/2];

DX = diff(x_new_nodes);
DY = diff(y_new_nodes);
%
D_ssd = bsxfun(@times,  [DX DY] , 1./sqrt(DX.^2+DY.^2));
Nssd = ( [0 -1; 1 0] *D_ssd')';
draw_ssd( fig , x_ls - Nssd*0.04 , Nssd, 0.035 );
% 


set(gca,'box','off')
color = get(fig,'Color');

set(gca,'XColor',color,'YColor',color,'TickDir','out')
set(gcf,'PaperPositionMode','auto');
print( fullfile( '../..','Figures/SFS_theory','WFS_ref_point' ) ,'-dpng')