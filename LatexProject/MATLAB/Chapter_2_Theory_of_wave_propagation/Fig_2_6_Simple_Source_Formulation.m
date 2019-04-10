clear;
close all
%
addpath(genpath('../Files'));

[pot_bb, pot_path] = read_epspath('general_contour.eps');
siz = norm(diff(pot_bb));
Le = siz/100;
radiator = meshpath(pot_path, siz/150);
radiator = translate_mesh(scale_mesh(radiator, 1/siz*5.25), [1.5 .5]);

Lx = 6.25;
Ly = 4;
field = create_slab([Lx, Ly], round([Lx Ly]*100));


[r_nodes, r_elem] = extract_core_mesh(radiator);
[f_nodes, f_elem] = extract_core_mesh(field);
f_elem(f_elem(:,1) == 32404,1) = 22404;

c = 343.1;
f0 = 1e3;
x0 = [0.4 2.5 0.];
k = 2*pi*f0/c;

[Ls, Ms, Lf, Mf] = helmholtz_bem_2d(r_nodes, r_elem, f_nodes, f_elem, k);

r_cent = centnorm(radiator);
ps = incident('line', x0, r_cent, [], k);

f_cent = centnorm(field);
pf_ana = incident('line', x0, f_cent, [], k);


I = speye(size(Ms));
qs_int = Ls \ (Ms - .5 * I) * ps;
qs_ext = Ls \ (Ms + .5 * I) * ps;

pf = Lf * (qs_ext - qs_int);
%%
ftsize = 9;
f = figure('Units','points','Position',[200,200,407,301]);
set(f,'defaulttextinterpreter','latex')

p1 = axes('Units','normalized','Position',[ 0.05 0.5 0.43 .6 ]);
hold on;
plot(x0(:,1), x0(:,2), 'k*');
plot_mesh(field, real(pf_ana));
shading flat;
plot_mesh(radiator);
axis equal tight;
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
caxis([-1 1] * 0.5e-1);

p2 = axes('Units','normalized','Position',[ 0.565 0.5 0.43 .6 ]);
hold on;
plot(x0(:,1), x0(:,2), 'k*');
plot_mesh(field, real(pf));
shading flat;
plot_mesh(radiator);
axis equal tight;
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
caxis([-1 1] * 0.5e-1);

p3 = axes('Units','normalized','Position',[ 0.31 -0.02  0.43 .6 ]);
hold on;
plot(x0(:,1), x0(:,2), 'k*');
plot_mesh(field, real(pf_ana - pf));
shading flat;
plot_mesh(radiator);
axis equal tight;
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
caxis([-1 1] * 0.5e-1);

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Basic_acoustics','simple_source_formulation' ) ,'-dpng')