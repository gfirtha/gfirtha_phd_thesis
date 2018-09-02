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
[ x_ssd, n_ssd, ds ] = create_equidist_SSD( [x_sec_contour y_sec_contour], Nc );

x_ref = (x_ssd(:,1)-mean(x_ssd(:,1)))*0.75 + mean(x_ssd(:,1));
y_ref = (x_ssd(:,2)-mean(x_ssd(:,2)))*0.75 + mean(x_ssd(:,2));
Diff_refX = x_ref(1:end-1)-x_ref(2:end); Diff_refY = y_ref(1:end-1)-y_ref(2:end);
D_ref = [ Diff_refX Diff_refY]; D_ref = bsxfun(@times, D_ref , 1./sqrt(Diff_refX.^2+Diff_refY.^2));
Nref_cont = ( [0 1; -1 0] *D_ref')';

reference_curve = [x_ref y_ref];
%reference_curve = reference_curve(1:end-1,:);

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

wb = waitbar(0,'Calculating');
for n = 1 : length(field_synth)
    waitbar(n/length(field_synth),wb);
    x_ =  cosd(phi)*( x_ssd(:,1)-xs ) + sind(phi)*(x_ssd(:,2)-ys);
    y_ = -sind(phi)*( x_ssd(:,1)-xs ) + cosd(phi)*(x_ssd(:,2)-ys);
    Tau_0 = sqrt( (x_ssd(:,1)-xr(n)).^2 + (x_ssd(:,2)-yr(n)).^2  )/c;
    R_0 = Tau_0*c;
    Delta = sqrt( ( x_-v*(t-Tau_0) ).^2 + ( y_.^2 ) * (1-M^2) );
    
    R = ( M*( x_-v*(t-Tau_0) ) + Delta )./(1-M^2);
    Tau = R / c;
    
    xse = cosd(phi)*v*(t-Tau_0-Tau)+xs;
    yse = sind(phi)*v*(t-Tau_0-Tau)+ys;
    K = [ x_ssd(:,1)-xse x_ssd(:,2)-yse  ];
    Kn = k*sum(K.*n_ssd,2)./Delta;
    Kn(Kn<0) = 0;
    
    
    TA = zeros(2,size(reference_curve,1)-1);
    dref = zeros(size(x_ssd,1),1);
    ref_win_x0  = zeros(size(x_ssd,1),1);
    for j = 1 : length(x_ssd)
     %   waitbar(j/length(x_ssd),wb);

        if Kn(j)~=0
        x0 = x_ssd(j,:)';
        V =  K(j,:)'./R(j);
        
        for m = 1 : size(reference_curve,1)-1
            x1 = reference_curve(m,:)';
            x2 = reference_curve(m+1,:)';
            A = [ V x2-x1 ];
            TA(:,m) = A\(x2-x0);
        end
        d = TA(1,:);
        a = TA(2,:);
        d( a<0 | a>1 ) = [];
        a( a<0 | a>1 ) = [];
        if ~isempty(d)
            dref(j) = min(d);
            ref_win_x0(j)  = 1;
        end         
        else end
    end
    x_ssd_refd = x_ssd(ref_win_x0 == 1,:);
    d_ref_valid = dref(ref_win_x0 == 1,:);
    K_refd = K(ref_win_x0==1,:)/k;

    tv = (1:length(dref))';
    dref_int = interp1(tv(ref_win_x0 == 1),dref(ref_win_x0 == 1,:),tv,'linear','extrap') ;
    int_t = 'nearest';
    dref_int = ( 0*interp1(x_ssd_refd(:,1),dref_int(ref_win_x0 == 1,:),x_ssd(:,1),int_t,'extrap')+...
         1*interp1(x_ssd_refd(:,2),dref_int(ref_win_x0 == 1,:),x_ssd(:,2),int_t,'extrap') );

     C = sqrt(2*pi*dref_int/(1i*k)).*sqrt(Delta./(dref_int+R));
    field_synth(n) = 1/(4*pi)*sum( 2*1i/(4*pi)*Kn.*C.*exp(1i*w0*(t-Tau-Tau_0))./(R_0.*Delta) ) *ds;
    
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
Kn = sum([ x_ssd(:,1)-xse x_ssd(:,2)-ys ].*n_ssd,2);
win0 = (Kn>0);
%%

ftsize = 13.75;
f = figure('Units','points','Position',[200,200,400,470]);
set(f,'defaulttextinterpreter','latex')
pos = [ 0.1   0.565  0.72 .45
    0.1   0.035 0.85 .5];

p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field_synth));
axis equal tight
caxis([-1 1] * 8e-2);
shading interp
hold on
line(x_sec_contour, y_sec_contour, z_sec_contour, 'Color', 'black','LineStyle',':','LineWidth',2);
line(x_ssd(win0,1), x_ssd(win0,2), 0*x_ssd(win0,1), 'Color', 'black','LineStyle','-','LineWidth',2);
plot( linspace(-3,3,100)*cosd(phi)+xs, linspace(-3,3,100)*sind(phi)+ys , ':w' , 'LineWidth',1)

%line(reference_curve_(:,1), reference_curve_(:,2), 0*reference_curve_(:,1), 'Color', 'white','LineStyle',':','LineWidth',1);
line(reference_curve(:,1), reference_curve(:,2), 0*reference_curve(:,1), 'Color', 'white','LineStyle',':','LineWidth',1);
%contour( x, y, real(field_ref),[0 0],'LineWidth',0.5,'LineColor',[1 1 1]*1);
%
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');

headWidth = 5;
headLength = 5;
LineLength = 0.75;
plot(xs,ys,'ok','MarkerSize',3,'MarkerFaceColor','black')
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
line(x_sec_contour, y_sec_contour, z_sec_contour, 'Color', 'black','LineStyle',':','LineWidth',2);
line(x_ssd(win0,1), x_ssd(win0,2), 0*x_ssd(win0,1), 'Color', 'black','LineStyle','-','LineWidth',2);
line(reference_curve(:,1), reference_curve(:,2), 0*reference_curve(:,1), 'Color', 'white','LineStyle',':','LineWidth',1);
plot( linspace(-3,3,100)*cosd(phi)+xs, linspace(-3,3,100)*sind(phi)+ys , ':w' , 'LineWidth',1)

xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
col = colorbar;
title(col,'[dB]'  , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
%
%
plot(xs,ys,'ok','MarkerSize',3,'MarkerFaceColor','black')
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[xs ys LineLength*cosd(phi) LineLength*sind(phi)]);
xlim([x(1),x(end)])
ylim([y(1),y(end)])

allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','25D_WFS_ms' ) ,'-dpng')
