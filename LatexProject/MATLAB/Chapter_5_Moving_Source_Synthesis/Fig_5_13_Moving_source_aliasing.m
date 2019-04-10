clear
close all
%%


c = 343.1;
v = c/4*2;
M = v / c;

fs = 50e3;
Lt =  0.3;
t = (-Lt:1/fs:Lt);
Nt = length(t);
w = fftshift( 2*pi*(-Nt/2:Nt/2-1)'/(Nt)*fs );

%Nf = 10;
%fir = fir1(Nf,[0.1,15000]/fs*2);
%s = circshift(filter(fir,1,fftshift([1; zeros(Nt-1,1)]),[]),-Nf/2);
L0 = 15;
s = [zeros(floor(Nt/2-L0/2),1);hann(L0);zeros(floor(Nt/2-L0/2),1)];

w0 = 2*pi*1000;
k = w0 / c;

dx = 0.10;

x_field = (-3:1e-2:3)'; 
y_field = (-1:1e-2:3)'; 
[X,Y] = meshgrid(x_field,y_field);

xs = [1.5,-2];
yref = 1.5;

s = circshift(s,-round(xs(1)/v*fs));
x0 = (-10:dx:10)';
y0 = zeros(size(x0));
x0 = [ x0 y0 ];


H = sqrt(1i*w/(2*pi*c));
d_wfs = real(ifft(H.*fft(s)));
aliased_field_freq_dom = zeros(size(X));
aliased_field_time_dom = zeros(size(X));
%%
t0 = (yref-xs(2))/c-xs(1)/v;
%t0 = 0.0125;
for n = 1 : size(x0,1)
    n
    r = sqrt((X-x0(n,1)).^2 + (Y-x0(n,2)).^2); 
    
    Delta = sqrt( (x0(n,1)-xs(1)-v*(t0-r/c) ).^2 + (0-    xs(2)).^2*(1-M^2) );
    R_dyn = ( M * ( x0(n,1)-xs(1) - v*(t0-r/c)  ) + Delta ) / ( (1-M^2) );
    
    aliased_field_time_dom = aliased_field_time_dom -...
            1/(4*pi)*sqrt(yref/(yref-xs(2)))*xs(2)*interp1(t,d_wfs,t0-(r+R_dyn)/c,'linear','extrap')./(r.*Delta.^(3/2))*dx;
    aliased_field_freq_dom = aliased_field_freq_dom - sqrt(1i*w0/(2*pi*c))*...
            1/(4*pi)*sqrt(yref/(yref-xs(2)))*xs(2)*exp(1i*w0*(t0-(r+R_dyn)/c))./(r.*Delta.^(3/2))*dx;
end

%%
ftsize = 9;
fig = figure('Units','points','Position',[200,200,407,144]);
set(fig,'defaulttextinterpreter','latex')

pos = [ 0.06   0.125  0.41 .9
        0.57   0.125  0.41 .9];

p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x_field,y_field,real(aliased_field_freq_dom));
shading interp
axis equal tight
caxis([-1,1]*5e-2)
caxis([-.05,.05])
hold on
%plot( [ -3+dx 3-dx ],  [ 0 0 ], 'k', 'Linewidth', 2 )
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
plot(x0(:,1),x0(:,2),'ok','MarkerSize',1.25,'MarkerFaceColor','black')
xlim([x_field(1), x_field(end)])

p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x_field,y_field,real(aliased_field_time_dom));
shading interp
axis equal tight
caxis([-1,1]*2e-2)
hold on
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
plot(x0(:,1),x0(:,2),'ok','MarkerSize',1.25,'MarkerFaceColor','black')
xlim([x_field(1), x_field(end)])

allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Moving_sources','Spatial_aliasing' ) ,'-dpng')
