clear
close all
%%
fs = 44.1e3;
Lt =  0.1;
t = (-Lt:1/fs:Lt);
Nt = length(t);
w = fftshift( 2*pi*(-Nt/2:Nt/2-1)'/(Nt)*fs );

Nf = 10;
fir = fir1(Nf,[10,10e3]/fs*2);
s = filter(fir,1,[1; zeros(Nt-1,1)],[]);

w0 = 2*pi*2e3;
c = 343.1;
k = w0 / c;

dx = 0.10;

x_field = (-3:1e-2:3)'; 
y_field = (-1:1e-2:3)'; 
[X,Y] = meshgrid(x_field,y_field);

xs = [0,-2];
yref = 1.5;

x0 = (-8:dx:8)';
y0 = zeros(size(x0));
x0 = [ x0 y0 ];

R = sqrt((x0(:,1)-xs(1)).^2+(x0(:,2)-xs(2)).^2);
H = sqrt(1i*w/(2*pi*c));
d_wfs = real(fftshift(ifft(H.*fft(s))));
Dxw = -1/(4*pi)*sqrt(8*pi/(1i*k))*sqrt(yref/(yref-xs(2)))*1i*k*xs(2)*exp(-1i*k*R)./R.^(3/2);
%%
aliased_field_freq_dom = zeros(size(X));
aliased_field_time_dom = zeros(size(X));

t0 = (yref-xs(2))/c+(Nf/2)/fs;
for n = 1 : size(x0,1)
    r = sqrt((X-x0(n,1)).^2 + (Y-x0(n,2)).^2); 
    aliased_field_time_dom = aliased_field_time_dom - 1/(4*pi)*sqrt(yref/(yref-xs(2)))*xs(2)*spline(t,d_wfs,t0-(R(n)+r)/c)./(r*R(n)^(3/2))*dx;
    aliased_field_freq_dom = aliased_field_freq_dom + Dxw(n)/(4*pi)*exp(-1i*k*r)./r*dx;
end

%%
ftsize = 14.3;
fig = figure('Units','points','Position',[200,200,650,230]);
set(fig,'defaulttextinterpreter','latex')   

pos = [ 0.06   0.125  0.41 .9
        0.57   0.125  0.41 .9];

p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x_field,y_field,real(aliased_field_freq_dom));
shading interp
axis equal tight
%caxis([-1,1]*5e-2)

caxis([-.05,.05])
hold on
%plot( [ -3+dx 3-dx ],  [ 0 0 ], 'k', 'Linewidth', 2 )
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
plot(x0(:,1),x0(:,2),'ok','MarkerSize',2,'MarkerFaceColor','black')
xlim([x_field(1), x_field(end)])

p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x_field,y_field,real(aliased_field_time_dom));
shading interp
axis equal tight
%caxis([-1,1]*5e-2)

caxis([-1,1]*1e-2)
hold on
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
plot(x0(:,1),x0(:,2),'ok','MarkerSize',2,'MarkerFaceColor','black')
xlim([x_field(1), x_field(end)])

    

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/SFS_theory','Spatial_alising' ) ,'-dpng')
