
clear;

%
c = 343.1;
f0 = 1e3;
xs = [0.4 2.5 0.];
k = 2*pi*f0/c;

[pot_bb, pot_path] = read_epspath('krumpli2.eps');
siz = norm(diff(pot_bb));
rad_int = meshpath(pot_path, siz/50);
rad_int = translate_mesh(scale_mesh(rad_int, 1/siz*5.25), [1.5 .5])
rad_ext = rad_int;
rad_ext.Elements = flip_elements(rad_ext.Elements);

Lx = 6.25;
Ly = 4;
field = create_slab([Lx, Ly], round([Lx Ly]*100));

%// extract the mesh description matrices
[f_nodes, f_elem] = extract_core_mesh(field);
%// convert slabs to 2D volume elements
f_elem(:,1) = 22404;

%% interior
[r_nodes, r_elem] = extract_core_mesh(rad_int);
[A, B, Lf_int, D] = helmholtz_bem_2d(r_nodes, r_elem, f_nodes, f_elem, k);

%% Generalized WFS

f_cent = centnorm(field);
P_ana = incident('line', xs, f_cent, [], k);

[C,n] = centnorm(rad_int);
ps_ana = incident('line', xs, C, [], k);
vx = bsxfun(@minus,C,xs);
vn = sum(vx.*n,2);
kn = k*(vn<=0)./sqrt(sum(vx.^2,2));


P_WFS = Lf_int * ( 2*1i*kn.*ps_ana ); % KH-approximation

%% Unified WFS 2D
ftsize = 14;

win = (vn<=0);
f = figure('Units','points','Position',[200,200,700,250]);

subplot(1,2,1);
plot_mesh(field, real( P_WFS ));
%set(gca, 'Units','normalized','Position',[ 0.05 0.1 0.42 .9 ]);
x = rad_int.Nodes(:,2);
y = rad_int.Nodes(:,3);
z = rad_int.Nodes(:,4);
line(x, y, z, 'Color', 'black','LineStyle',':','LineWidth',2);
line(rad_int.Nodes(win,2), rad_int.Nodes(win,3), rad_int.Nodes(win,4), 'Color', 'black','LineStyle','-','LineWidth',3);
caxis([-1 1] * 0.5e-1);
shading flat;
axis equal tight
set(gca,'box','on')
xlabel('$x \rightarrow \mathrm{[m]}$','Interpreter','latex','FontSize',ftsize)
ylabel('$z \rightarrow \mathrm{[m]}$','Interpreter','latex','FontSize',ftsize)
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');

subplot(1,2,2);
plot_mesh(field, real( P_ana - P_WFS ));
%set(gca, 'Units','normalized','Position',[ 0.05 0.1 0.42 .9 ]);
x = rad_int.Nodes(:,2);
y = rad_int.Nodes(:,3);
z = rad_int.Nodes(:,4);
line(x, y, z, 'Color', 'black','LineStyle',':','LineWidth',2);
line(rad_int.Nodes(win,2), rad_int.Nodes(win,3), rad_int.Nodes(win,4), 'Color', 'black','LineStyle','-','LineWidth',3);
caxis([-1 1] * 0.5e-1);
shading flat;
axis equal tight
set(gca,'box','on')
xlabel('$x \rightarrow \mathrm{[m]}$','Interpreter','latex','FontSize',ftsize)
ylabel('$z \rightarrow \mathrm{[m]}$','Interpreter','latex','FontSize',ftsize)
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');


set(gcf,'PaperPositionMode','auto');
%print -dpng KH_approx -r300