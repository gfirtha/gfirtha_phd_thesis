clear;
close all
addpath(genpath('../Files'));
%%
% Simulation parameters
Lx = 6.25;  % Field length
Ly = 4;     % Field height
dx = 1e-2;
fs = 50e3;

% Virtual source parameters
xs = [ 1.75 3 ];
Amp = 1.5;              % Virtual source amplitude             
source_exc = 'impulse';     % Excitation type: impulse / harmonic
loudspeaker_type = 'ideal'; % Loudspeaker directivity: ideal spatial LP ('ideal') / piston
ds = 0.1;

x = (0:dx:Lx);
y = (0:dx:Ly);
[X,Y] = meshgrid(x,y);
c = 343.1;
xrec0 = [3.25 1.5];
kv0 = -0.25;
%% Extract SSD mesh description matrices
[pot_bb, pot_path] = read_epspath('general_contour.eps');
siz = norm(diff(pot_bb));
rad_int = meshpath(pot_path, siz/150);
rad_int = translate_mesh(scale_mesh(rad_int, 1/siz*5.25), [1.5 .5]);
x_sec_contour = rad_int.Nodes(:,2);
y_sec_contour = rad_int.Nodes(:,3);
z_sec_contour = rad_int.Nodes(:,4);
%% Convert SSD into equidistant contour
L_ssd = sum(sqrt(diff(x_sec_contour).^2+diff(y_sec_contour).^2));
N_ssd = floor(L_ssd/ds);
[ x_ssd, N_ssd, ds ] = create_equidist_SSD( [x_sec_contour y_sec_contour], N_ssd);
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
K = bsxfun( @times, [x_ssd(:,1)-xs(1) x_ssd(:,2)-xs(2)] , 1./r0 );

kn = sum( K.*N_ssd, 2 );
kv = sum( K.*[N_ssd(:,2) -N_ssd(:,1)], 2 );
w0 = kn>=0;
x_ssd = x_ssd(w0,:);
N_ssd = N_ssd(w0,:);
K = K(w0,:);
kn =  kn(w0,:);
kv =  kv(w0,:);
r0 = r0(w0,:);
%% Find reference distance
[ dref,ref_win_x0 ] = find_ref_dist( x_ssd, reference_curve, K );
%%
x_ssd_refd = x_ssd(ref_win_x0 == 1,:);
d_ref_valid = dref(ref_win_x0 == 1,:);
K_refd = K(ref_win_x0==1,:);

t = (1:length(dref))';

dref = interp1(t(ref_win_x0 == 1),dref(ref_win_x0 == 1,:),t,'linear','extrap') ;
int_t = 'linear';
dref = ( 0.5*interp1(x_ssd_refd(:,1),dref(ref_win_x0 == 1,:),x_ssd(:,1),int_t,'extrap')+...
    0.5*interp1(x_ssd_refd(:,2),dref(ref_win_x0 == 1,:),x_ssd(:,2),int_t,'extrap') );
Dref = r0.*dref./(r0 + dref);
%%
L_t =  0.025;
t = (-L_t:1/fs:L_t);

s = zeros(length(t),1);
if strcmp(source_exc,'impulse')
    L0 = 5;
    s(round(end/2)+1:round(end/2)+L0) = hann(L0);
    s = s-mean(s);
elseif strcmp(source_exc,'harmonic')
    s = exp(1i*2*pi*2e3*t);
end
kn0 = 0;

[ field_synth ] = get_synthesized_field( t,s, Amp,xs,'off' ,loudspeaker_type, 'off', Dref,kn,r0, xrec0, x_ssd, N_ssd, X,Y,kv,kv0 );
[ field_synth_aa ] = get_synthesized_field( t,s, Amp,xs,'on' ,loudspeaker_type, 'off', Dref,kn,r0, xrec0, x_ssd, N_ssd, X,Y,kv,kv0 );
%
%%
ftsize = 9;
f = figure('Units','points','Position',[200,200,260,188]);
set(f,'defaulttextinterpreter','latex')
pos = [ 0.1   0.12  0.8 .85 ];


p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field_synth_aa));
axis equal tight
caxis([-1 1] * 2e-2);
shading interp
hold on
line(x_sec_contour, y_sec_contour, z_sec_contour, 'Color', 'black','LineStyle',':','LineWidth',1.25);
line(x_ssd(:,1), x_ssd(:,2), 0*x_ssd(:,1), 'Color', 'black','LineStyle','-','LineWidth',1.25);
plot(xs(1), xs(2), 'ok','MarkerSize',2, 'MarkerFaceColor','black')

[~,ind] =min(abs(kv-kv0));
plot(x_ssd(ind,1),x_ssd(ind,2),'.k','MarkerSize',12);
LineLength = 0.65;headWidth = 5*1.5;
headLength = 5*1.5;
ah = annotation('arrow','headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth,'LineWidth',1,'Color','Black','LineStyle','--');
set(ah,'parent',gca);
set(ah,'position',[x_ssd(ind,1),x_ssd(ind,2), LineLength*K(ind,1), LineLength*K(ind,2)]);    


reference_curve_ = [x_ref(~mask) y_ref(~mask)];
[~,i] = max(sqrt(sum(diff(reference_curve_,2).^2,2)));
reference_curve_ = [reference_curve_(i+2:end,:);reference_curve_(1:i+1,:)];
plot(x_ssd(:,1),x_ssd(:,2),'.k','MarkerSize',9);
xlabel( '$x$ [m]' , 'FontSize', ftsize );
ylabel( '$y$ [m]' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);

allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/SFS_theory','Steered_antialias_synth' ) ,'-dpng')