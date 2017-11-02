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
xs = [ 0.4 2.5 ];
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
[ x_ssd, N_ssd, ds ] = create_equidist_SSD( [x_sec_contour y_sec_contour], 1000 );

%% Create reference curve, with only the visible parts of the contour left
x_ref = (x_ssd(:,1)-mean(x_ssd(:,1)))*0.75 + mean(x_ssd(:,1));
y_ref = (x_ssd(:,2)-mean(x_ssd(:,2)))*0.75 + mean(x_ssd(:,2));
Diff_refX = x_ref(1:end-1)-x_ref(2:end); Diff_refY = y_ref(1:end-1)-y_ref(2:end);
D_ref = [ Diff_refX Diff_refY]; D_ref = bsxfun(@times, D_ref , 1./sqrt(Diff_refX.^2+Diff_refY.^2));
Nref_cont = ( [0 1; -1 0] *D_ref')'; 
kx = bsxfun ( @times , [ x_ref-xs(1) y_ref-xs(2) ] , 1./sqrt( (x_ref-xs(1)).^2 +  (y_ref-xs(2)).^2 ) );
mask =   ( sum( kx .*[Nref_cont; 0 0] , 2 ) >= 0 ) ;
reference_curve = [x_ref(mask) y_ref(mask)];
reference_curve = reference_curve(1:end-1,:);
%% Select the visible SSD elements
r0 = sqrt( (x_ssd(:,1)-xs(1)).^2+ ( x_ssd(:,2)-xs(2) ).^2 );
K = k*bsxfun( @times, [x_ssd(:,1)-xs(1) x_ssd(:,2)-xs(2)] , 1./r0 );
kn = sum( K.*N_ssd, 2 );
w0 = kn>=0;
x_ssd = x_ssd(w0,:);
N_ssd = N_ssd(w0,:);
K = K(w0,:);
kn =  kn(w0,:);
r0 = r0(w0,:);
%% Find reference distance
wb = waitbar(0,'Calculating');
TA = zeros(2,size(reference_curve,1)-1);
dref = zeros(size(x_ssd,1),1);
ref_win_x0  = zeros(size(x_ssd,1),1);
for n = 1 : length(x_ssd)
    waitbar(n/length(x_ssd),wb);

    x0 = x_ssd(n,:)';
    v =  K(n,:)'/k;
    
    for m = 1 : size(reference_curve,1)-1
        x1 = reference_curve(m,:)';
        x2 = reference_curve(m+1,:)';
        A = [ v x2-x1 ];
        TA(:,m) = A\(x2-x0);
    end
    d = TA(1,:);
    a = TA(2,:);
    d( a<0 | a>1 ) = [];
    a( a<0 | a>1 ) = [];
    if ~isempty(d)
       dref(n) = min(d); 
       ref_win_x0(n)  = 1;
    end
end
close(wb);
%%
x_ssd_refd = x_ssd(ref_win_x0 == 1,:);
d_ref_valid = dref(ref_win_x0 == 1,:);
K_refd = K(ref_win_x0==1,:)/k;

t = (1:length(dref))';

dref = interp1(t(ref_win_x0 == 1),dref(ref_win_x0 == 1,:),t,'linear','extrap') ;
int_t = 'nearest';
dref = ( 0*interp1(x_ssd_refd(:,1),dref(ref_win_x0 == 1,:),x_ssd(:,1),int_t,'extrap')+...
        1*interp1(x_ssd_refd(:,2),dref(ref_win_x0 == 1,:),x_ssd(:,2),int_t,'extrap') );
%int_t = 'linear';
%dref = ( 0.5*interp1(x_ssd_refd(:,1),dref(ref_win_x0 == 1,:),x_ssd(:,1),int_t,'extrap')+...
%         0.5*interp1(x_ssd_refd(:,2),dref(ref_win_x0 == 1,:),x_ssd(:,2),int_t,'extrap') );
%
Dref = r0.*dref./(r0 + dref);

Dx0 = 1i/(4*pi)* sqrt(8*pi/(1i*k)).*sqrt(Dref).*kn.*exp(-1i*k*r0)./r0;
%%
field_synth = zeros(size(X));

wb = waitbar(0,'Calculating');
for n = 1 : length(Dx0)
   
    waitbar(n/length(Dx0),wb);
    R = sqrt( (X-x_ssd(n,1)).^2 + (Y-x_ssd(n,2)).^2 );
    field_synth = field_synth + 1/(4*pi)*Dx0(n)*exp(-1i*k*R)./R*ds;
   
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
line(x_ssd(:,1), x_ssd(:,2), 0*x_ssd(:,1), 'Color', 'black','LineStyle','-','LineWidth',2);

reference_curve_ = [x_ref(~mask) y_ref(~mask)];
[~,i] = max(sqrt(sum(diff(reference_curve_,2).^2,2)));
reference_curve_ = [reference_curve_(i+2:end,:);reference_curve_(1:i+1,:)];
line(reference_curve_(:,1), reference_curve_(:,2), 0*reference_curve_(:,1), 'Color', 'white','LineStyle',':','LineWidth',1);
line(reference_curve(:,1), reference_curve(:,2), 0*reference_curve(:,1), 'Color', 'white','LineStyle','--','LineWidth',1);

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
line(x_ssd(:,1), x_ssd(:,2), 0*x_ssd(:,1), 'Color', 'black','LineStyle','-','LineWidth',2);


line(reference_curve_(:,1), reference_curve_(:,2), 0*reference_curve_(:,1), 'Color', 'white','LineStyle',':','LineWidth',1);
line(reference_curve(:,1), reference_curve(:,2), 0*reference_curve(:,1), 'Color', 'white','LineStyle','--','LineWidth',1);   
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
c = colorbar;
title(c,'[dB]' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
%
%
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( fullfile( '../..','Figures/SFS_theory','25D_WFS_general' ) ,'-dpng')