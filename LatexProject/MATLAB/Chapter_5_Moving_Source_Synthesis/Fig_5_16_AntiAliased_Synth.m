clear
close all
addpath(genpath('../Files'));

c = 343.1;
v = c*3/4;


N_ssd = 125;                                 % Number of loudspeakers
fi = (0:2*pi/N_ssd:2*pi-2*pi/N_ssd)';
Rssd = 2;
x0 = Rssd*[ cos(fi) sin(fi) ];
dx0 = mean(sqrt(sum(diff(x0,1).^2,2)));

Rref = 0.1;
xref = Rref*[ cos(fi) sin(fi) ];

n0 = -[ cos(fi) sin(fi) ];
v0 =  [ sin(fi) -cos(fi)];
%
dx = .5e-2;
x_field = (-3.5:dx:2.5);
y_field = (-2.5:dx:3.5 );
[X,Y] = meshgrid(x_field,y_field);

%
x_a = [   -3   -3  -3  -2.5   -1.5    0  2.5;    %x_coordinates
    -2.5  -1   0   1.5    2.5    3  3 ];  %y_coordinates
[ p,xp,yp ] = make_path( x_a(1,:), x_a(2,:), 150 );

Tsim = p(end)/v;
%Tsim = 0.02;
fs = 50e3;
t = (0:1/fs:Tsim-1/fs)';
Nt = length(t);
w = 2*pi*fftshift( (-Nt/2:Nt/2-1)/(Nt)*fs );
L0 = 9;
n1 = 50; n2 = 3*n1;
s0 = [ zeros(n1,1);hann(L0);zeros(n2,1)];
s = repmat(s0,ceil(length(t)/length(s0)),1)';
s = s(1:length(t));
imp_pos = n1+ceil(L0/2)+(0:ceil(length(t)/length(s0)))*(n1+n2+L0);

xs = get_trajectory( p,xp,yp, t, v );
%% Calculate driving functions
Nt = length(t);
w = 2*pi*fftshift( (-Nt/2:Nt/2-1)/(Nt)*fs );
d_wfs = zeros(length(t),length(x0));
s_wfs = ifft( sqrt(1i*w/(c*2*pi)).*fft(s) );
%
Tau = get_initial_position( v,c, x_a, x0 );
wb = waitbar(0,'Calculating driving functions');

wc = zeros(length(t),length(x0));
for n = 1 : length(t)
    waitbar(n/length(t),wb);
    xs_t = interp1( t,xs(:,1), t(n)-Tau, 'linear','extrap' );
    ys_t = interp1( t,xs(:,2), t(n)-Tau, 'linear','extrap' );
    vx = (xs_t - interp1( t,xs(:,1), t(n)-1/fs-Tau, 'linear','extrap' ))*fs;
    vy = (ys_t - interp1( t,xs(:,2), t(n)-1/fs-Tau, 'linear','extrap' ))*fs;
    Dvx = x0-[xs_t ys_t];
    R = sqrt( sum( Dvx.^2,2) );
    Vv = 1/c*sum([vx vy].*(Dvx),2);
    Delta = R - Vv;
    D = xs_t.*x0(:,2) - ys_t.*x0(:,1);
    xref = ( D.*Dvx(:,2) - sign(Dvx(:,2)).*Dvx(:,1).*sqrt(Rref^2.*R.^2-D.^2 ))./R.^2;
    yref = (-D.*Dvx(:,1) - abs(Dvx(:,2)).*sqrt(Rref^2.*R.^2-D.^2 ))./R.^2;
    dref = sqrt( sum( (x0-[xref yref]).^2, 2)  );
    Kn = sum(n0.*Dvx ,2);
    Kn =  Kn.*(Kn>0);
    d_wfs(n,:) = sqrt(dref./(dref+R)).*Kn.*interp1(t,s_wfs,t(n)-Tau,'linear')./Delta.^(3/2)*dx0;
    kt = sum(Dvx./R.*v0,2);
    wc(n,:) = pi/dx0*c./abs(kt);
    Tau = Tau - 1/fs*Vv./Delta;
end
close(wb);
d_wfs(isnan(d_wfs)) = 0;
%% Antialising filtering
wlen = 20;
d_wfs_aa  = anti_aliasing_td( d_wfs, t, t,wc, wlen, wlen);
%%
Rfield_full = reshape( sqrt((bsxfun( @minus, X(:), x0(:,1)' )).^2+ (bsxfun( @minus, Y(:), x0(:,2)' )).^2),...
    length(y_field), length(x_field) , length(x0));

t0 = 0.02;
field_synth    = zeros(size(X));
field_synth_aa = zeros(size(X));
for n = 1 : length(x0)
    Rfield_synth = squeeze( Rfield_full(:,:,n) );
    field_synth    = field_synth    + 1/(4*pi)*interp1( t, d_wfs(:,n)   , t0-Rfield_synth/c,'nearest','extrap'   )./Rfield_synth;
    field_synth_aa = field_synth_aa + 1/(4*pi)*interp1( t, d_wfs_aa(:,n), t0-Rfield_synth/c,'nearest','extrap'   )./Rfield_synth;
end

%%
[~,i] = min(abs(t-t0));
V = (xs(i+1,:)-xs(i-1,:))*fs/2;

ftsize = 20;
f = figure('Units','points','Position',[200,100,850,370]);
set(f,'defaulttextinterpreter','latex')

fig_pos = [ 0.0   0.135  0.52 .85;
    0.5   0.135 0.52 .85];

sp1 = axes('Units','normalized','Position',fig_pos(1,:));
p1 = pcolor(sp1, x_field,y_field,real(field_synth));
shading interp
axis equal tight
caxis([-1 1] * 1.5e-2);
hold on
plot(sp1, xs(:,1),xs(:,2),':w','LineWidth',2)
xlim([x_field(1),x_field(end)])
ylim([y_field(1),y_field(end)])
%draw_ssd( sp1, x0(1:1:end,:), n0(1:1:end,:), 0.03 )

plot( x0(:,1), x0(:,2), 'ok','MarkerSize',3, 'MarkerFaceColor','black')
%plot(xs(1), xs(2), 'ok','MarkerSize',3, 'MarkerFaceColor','black')

pos1 = plot( sp1, xs(i,1),xs(i,2),'ok','MarkerFaceColor','black');
plot( sp1, xs(imp_pos(1:5),1),xs(imp_pos(1:5),2),'ok','MarkerFaceColor','white');
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);


headWidth = 5;
headLength = 5;
LineLength = 0.75;
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[xs(i,1) xs(i,2) LineLength*V(1)/v LineLength*V(2)/v]);


sp2 = axes('Units','normalized','Position',fig_pos(2,:));
p2 = pcolor(sp2, x_field,y_field,real(field_synth_aa));
shading interp
axis equal tight
caxis([-1 1] * 1.5e-2);
hold on
plot(sp2, xs(:,1),xs(:,2),':w','LineWidth',2)
pos2 = plot( sp2, xs(i,1),xs(i,2),'ok','MarkerFaceColor','black');
plot( sp2, xs(imp_pos(1:5),1),xs(imp_pos(1:5),2),'ok','MarkerFaceColor','white');
xlim([x_field(1),x_field(end)])
ylim([y_field(1),y_field(end)])
%draw_ssd( sp2, x0(1:1:end,:), n0(1:1:end,:), 0.03 )
plot( x0(:,1), x0(:,2), 'ok','MarkerSize',3, 'MarkerFaceColor','black')
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[xs(i,1) xs(i,2) LineLength*V(1)/v LineLength*V(2)/v]);

%%

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','antialiased_synth_moving_source' ) ,'-dpng')
