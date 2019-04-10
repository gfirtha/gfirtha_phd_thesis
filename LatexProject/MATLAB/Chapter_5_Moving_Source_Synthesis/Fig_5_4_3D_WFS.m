clear;
close all
addpath(genpath('../Files'));
Lx = 6.25;
Ly = 4;
dx = 2e-2;
x = (0:dx:Lx);
y = (0:dx:Ly);
[X,Y] = meshgrid(x,y);

c = 343.1;
f0 = 1.5e3;
w0 = 2*pi*f0;
k = w0/c;

% Extract SSD mesh description matrices
[pot_bb, pot_path] = read_epspath('general_contour.eps');
siz = norm(diff(pot_bb));
rad_int = meshpath(pot_path, siz/150);
rad_int = translate_mesh(scale_mesh(rad_int, 1/siz*5.25), [1.5 .5]);
x_sec_contour = rad_int.Nodes(:,2);
y_sec_contour = rad_int.Nodes(:,3);
z_sec_contour = rad_int.Nodes(:,4);
%% Convert SSD into equidistant contour
Nc = 800;
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


%%
field_synth = zeros(size(X));
field_synth = field_synth(:);
xr = X(:);
yr = Y(:);
zr = 0*Y(:);

xs = 1.85; ys = 3.2; zs = 0;

v = c/2;
M = v/c;
t = 0;

phi = 35;

%tic
p = ProgressBar(length(field_synth));
parfor n = 1:length(field_synth)
%wb = waitbar(0,'Calculating');
%for n = 1 : length(field_synth)
%   waitbar(n/length(field_synth),wb);
   x_ =  cosd(phi)*( X_SSD(:,1)-xs ) + sind(phi)*(X_SSD(:,2)-ys);
   y_ = -sind(phi)*( X_SSD(:,1)-xs ) + cosd(phi)*(X_SSD(:,2)-ys);
     p.progress ;  
   Tau_0 = sqrt( (X_SSD(:,1)-xr(n)).^2 + (X_SSD(:,2)-yr(n)).^2 + (X_SSD(:,3)-zr(n)).^2 )/c;
   R_0 = Tau_0*c;
   Delta = sqrt( ( x_-v*(t-Tau_0) ).^2 + ( y_.^2 + ( X_SSD(:,3)-zs ).^2 ) * (1-M^2) );
   
   R = ( M*( x_-v*(t-Tau_0) ) + Delta )./(1-M^2);
   Tau = R / c;
   
   xse = cosd(phi)*v*(t-Tau_0-Tau)+xs;
   yse = sind(phi)*v*(t-Tau_0-Tau)+ys;
   K = [ X_SSD(:,1)-xse X_SSD(:,2)-yse X_SSD(:,3)-zs ];
   Kn = k*sum(K.*N_SSD,2)./Delta;
   Kn(Kn<0) = 0;
   field_synth(n) = 1/(4*pi)*sum( 2*1i/(4*pi)*Kn.*exp(1i*w0*(t-Tau-Tau_0))./(R_0.*Delta) ) *ds*dz;
   
end
close(wb);
%toc
field_synth = reshape(field_synth,length(y),length(x));

%%
X_ = cosd(phi)*(X-xs)+sind(phi)*(Y-ys);
Y_ = -sind(phi)*(X-xs)+cosd(phi)*(Y-ys);
[X,Y] = meshgrid(x,y);

Delta_0 = sqrt( ( X_ - v*t ).^2 + ( Y_ ).^2*(1-M^2));
R_0 = (M*( X_ - v*t ) + Delta_0)/(1-M^2);

field_ref = 1/(4*pi)*exp(1i*w0*(t-R_0/c))./Delta_0;

%%

Delta_0 = sqrt( (x_ssd(:,1)-xs-v*t).^2 +  (x_ssd(:,2)-ys).^2*(1-M^2) );
R_0 = (M*(x_ssd(:,1)-xs-v*t) + Delta_0 )./(1-M^2);
tau_0 = R_0/c;
xse = v*(t-tau_0)+xs;
Kn = sum([ x_ssd(:,1)-xse x_ssd(:,2)-ys ].*N_ssd,2);
win0 = (Kn>0);
%%
ftsize = 9;
f = figure('Units','points','Position',[200,200,260,306]);
set(f,'defaulttextinterpreter','latex')
pos = [ 0.1   0.565  0.72 .45
        0.1   0.035 0.89 .5];
    
    
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field_synth));
axis equal tight
caxis([-1 1] * 8e-2);
shading interp
hold on
line(x_sec_contour, y_sec_contour, z_sec_contour, 'Color', 'black','LineStyle',':','LineWidth',1.25);
line(x_ssd(win0,1), x_ssd(win0,2), 0*x_ssd(win0,1), 'Color', 'black','LineStyle','-','LineWidth',1.25);
plot( linspace(-3,3,100)*cosd(phi)+xs, linspace(-3,3,100)*sind(phi)+ys , ':w' , 'LineWidth',1)
%contour( x, y, real(field_ref),[0 0],'LineWidth',0.5,'LineColor',[1 1 1]*1);
% 
xlabel( '$x$ [m]' , 'FontSize', ftsize );
ylabel( '$y$ [m]' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');

headWidth = 5;
headLength = 5;
LineLength = 0.75;
plot(xs,ys,'ok','MarkerSize',2,'MarkerFaceColor','black')
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[xs ys LineLength*cosd(phi) LineLength*sind(phi)]);
xlim([x(1),x(end)])
ylim([y(1),y(end)])
%
%
%
p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x,y, 20*log10(abs(field_ref - field_synth)));
axis equal tight
caxis([-50 ,-20]);
shading interp
hold on
line(x_sec_contour, y_sec_contour, z_sec_contour, 'Color', 'black','LineStyle',':','LineWidth',1.25);
line(x_ssd(win0,1), x_ssd(win0,2), 0*x_ssd(win0,1), 'Color', 'black','LineStyle','-','LineWidth',1.25);
plot( linspace(-3,3,100)*cosd(phi)+xs, linspace(-3,3,100)*sind(phi)+ys , ':w' , 'LineWidth',1)

xlabel( '$x$ [m]' , 'FontSize', ftsize );
ylabel( '$y$ [m]' , 'FontSize', ftsize );
col = colorbar;
title(col,'[dB]'  , 'FontSize', ftsize-1);
set(gca,'FontName','Times New Roman');
%
%
plot(xs,ys,'ok','MarkerSize',2,'MarkerFaceColor','black')
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[xs ys LineLength*cosd(phi) LineLength*sind(phi)]);
xlim([x(1),x(end)])
ylim([y(1),y(end)])

allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','3D_WFS2' ) ,'-dpng')
