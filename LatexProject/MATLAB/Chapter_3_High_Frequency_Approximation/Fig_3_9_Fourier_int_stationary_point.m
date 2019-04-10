clear
close all

c = 343.1;
w = 2*pi*1e3;
k = w/c;
lambda = 2*pi/k;

L = 2;
dx = 0.01;
x = (-L:dx:L)';
y = (-1:dx:L-0.5)';

xs = [0,-0.75];
[X,Y] = meshgrid(x,y);

kx = 0.5*k;

nkx = (x-xs(1))./sqrt( ( x-xs(1) ).^2 + xs(2).^2 );
nky = -xs(2)./sqrt( (x-xs(1)).^2 + xs(2).^2 );
[~,ind] = min( (nkx-kx/k).^2 );

x0 = [x(ind) 0];
phi0 = -1i*(k*norm( x0-xs ) + pi/2);

r = sqrt((X-xs(1)).^2+(Y-xs(2)).^2);
Pps = 1/(4*pi)*exp(-1i*k*r-phi0)./r;

phi0pw = 1i*(kx*x0(1) + pi/2);
Ppw = exp(1i*(kx*X+sqrt(k^2-kx^2)*Y)-phi0pw);

Nx = 5e3;
kx0 = 2*pi*(-Nx/2:Nx/2-1)'/(Nx*dx);
Pps_spec = -1i/4*besselh(0,2,-1i*sqrt(kx0.^2-k^2)*abs(xs(2)));
%%
fig = figure('Units','points','Position',[200,200,407,220]);
set(fig,'defaulttextinterpreter','latex')
ftsize = 9;

pos = [0.01 0.34 .53  1 ;
       0.01 -0.165 .53 1 ;
       0.66 0.635  .325 .325;
       0.66 0.125  .325 .325 ];

p1 = axes('Units','normalized','Position',pos(1,:));

ind2 = find((-0.25<y).*(y<0.75));
y1 = y(ind2);
p1_ = pcolor( x,y1, real(Ppw(ind2,:)) );
set(p1_, 'EdgeColor', 'white');

shading interp
axis equal tight
caxis([-1,1]*2e0)
hold on
contour( x, y1, real(Ppw(ind2,:)),1, '-k');
set(gca,'box','off')
color = get(fig,'Color');
set(gca,'XColor',color,'YColor',color,'TickDir','out')
ol = 15*dx;
headWidth = 5;
headLength = 5;
LineLength = 0.7;
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[x(1)-ol, 0, x(end)-x(1)+ol*2, 0]);    

ah2 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah2,'parent',gca);
set(ah2,'position',[-1.5, y1(1)-ol, 0, y1(end)-y1(1)+ol*2]);  

xlim([x(1)-ol,x(end)+ol]);
ylim([y1(1)-ol,y1(end)+ol]);


q = 15;

for i = rem(ind,q):q:length(x)
        ah = annotation('arrow',...
            'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[x(i) 0 LineLength*nkx(ind) LineLength*nky(ind)]);
        if ( i == ind )
            set(ah,'color','white');
        end
end

plot(x0(1),x0(2),'.k','MarkerSize',10);

p2 = axes('Units','normalized','Position',pos(2,:));

p2_ = pcolor( x,y, real(Pps) );
set(p2_, 'EdgeColor', 'white');

shading interp
axis equal tight
caxis([-1,1]*1e-1)
hold on
pcolor_ax =(gca);
cRange = caxis;
[C,~]= contour( x, y, real(Pps), [ 0 0 ], '-k');

r_t = sqrt( (X-x0(1)).^2 + (Y-x0(2)).^2 );
Ppw( r_t> 1 ) = nan;
[C2,~]= contour( x, y, real(Ppw), 1,'--k');

ol = 15*dx;
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[x(1)-ol, 0, x(end)-x(1)+ol*2, 0]);    

ah2 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah2,'parent',gca);
set(ah2,'position',[-1.5, y(1)-ol, 0, y(end)-y(1)+ol*2]);  

xlim([x(1)-ol,x(end)+ol]);
ylim([y(1)-ol,y(end)+ol]);

q = 15;
for i = rem(ind,q):q:length(x)
        ah = annotation('arrow',...
            'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[x(i) 0 LineLength*nkx(i) LineLength*nky(i)]);
        if ( i == ind )
            set(ah,'color','white');
        end
end
line([x(1),x(end)],[0,0],'Color','black','LineStyle','-','LineWidth',2)
plot(xs(1),xs(2),'.k','MarkerSize',15)
plot(x0(1),x0(2),'.k','MarkerSize',15)
set(gca,'box','off')
color = get(fig,'Color');
set(gca,'XColor',color,'YColor',color,'TickDir','out');

p3 = axes('Units','normalized','Position',pos(3,:));
p3_ = plot( x/lambda, nkx,'LineWidth',1);
xl = [x(1) x(end)]/lambda;
yl = ylim;
xlabel( '$x/\lambda$' , 'FontSize', ftsize );
ylabel( '$\hat{k}_{x}^P(x,0,0)$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
grid on
hold on
plot(x(ind)/lambda,nkx(ind),'ko','MarkerSize',3,'MarkerFaceColor','black')
line([xl(1),x(ind)/lambda],[nkx(ind),nkx(ind)],'Color','black','LineStyle','--');
line([x(ind),x(ind)]/lambda,[yl(1),nkx(ind)],'Color','black','LineStyle','--');
xlim(xl);

p4 = axes('Units','normalized','Position',pos(4,:));
plot(kx0/k,abs(Pps_spec),'LineWidth',1);
xlabel( '$k_x/k$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
ylabel( '$|\tilde{P}(k_x,0,0)|$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
grid on
hold on
Pkx0 = -1i/4*besselh(0,2,-1i*sqrt(kx.^2-k^2)*abs(xs(2)));
plot(kx/k,abs( Pkx0 ), 'ko','MarkerSize',3,'MarkerFaceColor','black')
line([kx/k,kx/k],[0,abs(Pkx0)],'Color','black','LineStyle','--');

xlim(1.2*[-1,1])
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/High_freq_approximations','fourier_stat_point' ) ,'-dpng')