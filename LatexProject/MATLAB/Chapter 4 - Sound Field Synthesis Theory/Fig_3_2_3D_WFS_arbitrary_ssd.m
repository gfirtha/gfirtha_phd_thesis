clear;
close all
addpath(genpath('../Files'));
Lx = 6.25;
Ly = 4;
dx = 1e-2;
x = (0:dx:Lx);
y = (0:dx:Ly);
[X,Y] = meshgrid(x,y);

c = 343.1;
f0 = 1.5e3;
k = 2*pi *f0/c;
xs = [ 0.4 2.5 0 ];
field_ref = 1/(4*pi)*exp( -1i*k* sqrt( (X-xs(1)).^2 + (Y-xs(2)).^2) )./sqrt( (X-xs(1)).^2 + (Y-xs(2)).^2);


% Extract SSD mesh description matrices
[pot_bb, pot_path] = read_epspath('general_contour.eps');
siz = norm(diff(pot_bb));
rad_int = meshpath(pot_path, siz/150);
rad_int = translate_mesh(scale_mesh(rad_int, 1/siz*5.25), [1.5 .5]);
x_sec_contour = rad_int.Nodes(:,2);
y_sec_contour = rad_int.Nodes(:,3);
z_sec_contour = rad_int.Nodes(:,4);
%% Convert SSD into equidistant contour
Nc = 700;
Nz = 401;
[ x_ssd, N_ssd, ds ] = create_equidist_SSD( [x_sec_contour y_sec_contour], Nc );
z_ssd = linspace(-4,4,Nz)';
dz = z_ssd(2)-z_ssd(1);
win = tukeywin(length(z_ssd),0.25);
z_ssd = repmat(z_ssd,1,Nc)';
z_ssd=z_ssd(:);
win = repmat(win,1,Nc)';
win=win(:);
X_SSD = [repmat(x_ssd,Nz,1) z_ssd];
N_SSD = [repmat(N_ssd,Nz,1) zeros(size(z_ssd,1),1)];

%% Select the visible SSD elements
r0 = sqrt( sum ( bsxfun( @minus, X_SSD , xs ) .^2 , 2) );

K = k*bsxfun( @times, bsxfun( @minus, X_SSD , xs ) , 1./r0 );

kn = sum( K.*N_SSD, 2 );
w0 = kn>=0;
X_SSD = X_SSD(w0,:);
N_SSD = N_SSD(w0,:);
K = K(w0,:);
kn =  kn(w0,:);
r0 = r0(w0,:);
win = win(w0);

Dx0 = 1/(4*pi)*2*1i.*kn.*exp(-1i*k*r0)./r0.*win;
%%
field_synth = zeros(size(X));

wb = waitbar(0,'Calculating');
for n = 1 : length(Dx0)
   
    waitbar(n/length(Dx0),wb);
    R = sqrt( (X-X_SSD(n,1)).^2 + (Y-X_SSD(n,2)).^2 + X_SSD(n,3).^2 );
    field_synth = field_synth + 1/(4*pi)*Dx0(n)*exp(-1i*k*R)./R*ds*dz;
   
end
close(wb)
%%
ftsize = 13.75;
f = figure('Units','points','Position',[200,200,400,470]);
set(f,'defaulttextinterpreter','latex')
pos = [ 0.1   0.565  0.72 .45
        0.1   0.035 0.85 .5];
    
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field_synth));
axis equal tight
caxis([-1 1] * 5e-2);
shading interp
hold on
line(x_sec_contour, y_sec_contour, z_sec_contour, 'Color', 'black','LineStyle',':','LineWidth',2);
line(x_ssd(w0(1:Nc),1), x_ssd(w0(1:Nc),2), 0*x_ssd(w0(1:Nc),1), 'Color', 'black','LineStyle','-','LineWidth',2);

xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');

%
%
%
p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x,y,20*log10(abs(field_ref - field_synth)));
axis equal tight
caxis([-52 ,-20]);
shading interp
hold on
line(x_sec_contour, y_sec_contour, z_sec_contour, 'Color', 'black','LineStyle',':','LineWidth',2);
line(x_ssd(w0(1:Nc),1), x_ssd(w0(1:Nc),2), 0*x_ssd(w0(1:Nc),1), 'Color', 'black','LineStyle','-','LineWidth',2);

xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
c = colorbar;
title(c,'[dB]'  , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
%
%
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( fullfile( '../..','Figures/SFS_theory','3D_WFS_general' ) ,'-dpng')