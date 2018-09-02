clear;
close all
%
addpath(genpath('../Files'));

[pot_bb, pot_path] = read_epspath('general_contour.eps');
siz = norm(diff(pot_bb));
radiator = meshpath(pot_path, siz/500);
radiator = translate_mesh(scale_mesh(radiator, 1/siz*5.25), [1.5 .5]);

Lx = 6.25;
Ly = 4;
field = create_slab([Lx, Ly], round([Lx Ly]*100));


[r_nodes, r_elem] = extract_core_mesh(radiator);
[f_nodes, f_elem] = extract_core_mesh(field);
f_elem(f_elem(:,1) == 32404,1) = 22404;

c = 343.1;
f0 = 1e3;
xs = [ 0.4 2.5 0.];
k = 2*pi*f0/c;

[Ls, Ms, Lf, Mf] = helmholtz_bem_2d(r_nodes, r_elem, f_nodes, f_elem, k);

r_cent = centnorm(radiator);
ps = incident('line', xs, r_cent, [], k);

f_cent = centnorm(field);
pf_ana = incident('line', xs, f_cent, [], k);

[x0,n] = centnorm(radiator);
ps_ana = incident('line', xs, x0, [], k);
vx = bsxfun(@minus,x0,xs);
vn = sum(vx.*n,2);
win = (vn<=0);

I = speye(size(Ms));
qs_int = Ls \ (Ms - .5 * I) * ps;
qs_ext = Ls \ (Ms + .5 * I) * ps;

pf = Lf * ( -2* qs_int.*win);
%%
ftsize = 11;
f = figure('Units','points','Position',[200,200,500,180]);
set(f,'defaulttextinterpreter','latex')

p1 = axes('Units','normalized','Position',[ 0.015 0.2 0.5 .75]);
hold on;
plot_mesh(field, real(pf));
shading flat;
plot(x0(:,1), x0(:,2), ':k','LineWidth',2);
plot(x0(win,1), x0(win,2), 'k','LineWidth',2);
axis equal tight;
xlabel( '$x \rightarrow$ [m]', 'FontSize', ftsize );
ylabel( '$y \rightarrow$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
caxis([-1 1] * 0.5e-1);

p2 = axes('Units','normalized','Position',[ 0.525 0.2 0.5 .75 ]);
hold on;
plot(xs(:,1), xs(:,2), 'k*');
plot_mesh(field, real(-pf+pf_ana));
shading flat;
plot_mesh(radiator);
axis equal tight;
xlabel( '$x \rightarrow$ [m]', 'FontSize', ftsize );
ylabel( '$y \rightarrow$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
caxis([-1 1] * 0.5e-1);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/High_freq_approximations','KH_approx.png' ) ,'-dpng')