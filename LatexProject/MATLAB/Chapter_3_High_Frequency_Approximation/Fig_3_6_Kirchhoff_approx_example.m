clear;
close all
%
addpath(genpath('../Files'));

[pot_bb, pot_path] = read_epspath('general_contour.eps');
siz = norm(diff(pot_bb));
radiator = meshpath(pot_path, siz/200);
radiator = translate_mesh(scale_mesh(radiator, 1/siz*5.25), [1.5 .5]);
%%
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
tic
[Ls, Ms, Lf, Mf] = helmholtz_bem_2d(r_nodes, r_elem, f_nodes, f_elem, k);
toc
[r_cent, r_norm] = centnorm(radiator);
[ps,qs] = incident('line', xs, r_cent, r_norm, k);

f_cent = centnorm(field);
pf_ana = incident('line', xs, f_cent, [], k);

[x0,n] = centnorm(radiator);
ps_ana = incident('line', xs, x0, [], k);
vx = bsxfun(@minus,x0,xs);
vn = sum(vx.*n,2);
win = (vn<=0);

pf = Lf * ( 2 * qs.*win);
%%
f = figure('Units','points','Position',[200,200,500,370]);
set(f,'defaulttextinterpreter','latex')

p1 = axes('Units','normalized','Position',[ 0.05 0.5 0.43 .6 ]);
hold on;
plot_mesh(field, real(pf));
shading flat;
plot( x0(:,1),   x0(:,2), ':k','LineWidth',2);
plot( x0(win,1), x0(win,2), 'k','LineWidth',2);
plot( xs(:,1),   xs(:,2), 'ok','MarkerSize',2, 'MarkerFaceColor','black')
axis equal tight;
xlabel( '$x \rightarrow$ [m]', 'FontSize', ftsize );
ylabel( '$y \rightarrow$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
caxis([-1 1] * 0.5e-1);

p2 = axes('Units','normalized','Position',[ 0.565 0.5 0.43 .6 ]);
hold on;
plot(xs(:,1), xs(:,2), 'k*');
plot_mesh(field, real(-pf+pf_ana));
shading flat;
plot_mesh(radiator);
plot( xs(:,1), xs(:,2), 'ok','MarkerSize',2, 'MarkerFaceColor','black')
axis equal tight;
xlabel( '$x \rightarrow$ [m]', 'FontSize', ftsize );
ylabel( '$y \rightarrow$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
caxis([-1 1] * 0.5e-1);

p3 = axes('Units','normalized','Position',[ 0.25 -0.02  0.55 .6 ]);
hold on;
plot(xs(:,1), xs(:,2), 'k*');
plot_mesh(field, 20*log10(abs(-pf+pf_ana)));
shading flat;
plot_mesh(radiator);
plot( xs(:,1), xs(:,2), 'ok','MarkerSize',2, 'MarkerFaceColor','black')
axis equal tight;
xlabel( '$x \rightarrow$ [m]', 'FontSize', ftsize );
ylabel( '$y \rightarrow$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
caxis([-50,-20]);
c = colorbar;
title(c,'[dB]' , 'FontSize', ftsize-1);
set(gca,'FontName','Times New Roman');

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/High_freq_approximations','KH_approx.png' ) ,'-dpng')