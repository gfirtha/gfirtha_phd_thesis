
clear 
close

dx = 0.01;

L = 2;
x = (-L:dx:L)';
y = (-0:dx:L)';
z = (-0:dx:L)';

[Xxy,Yxy] = meshgrid(x,y);
[Xxz,Zxz] = meshgrid(x,z);
[Yyz,Zyz] = meshgrid(y,z);


c = 343.1;
w = 2*pi*1e3;
k = w/c;
fi = 60;
theta = 90;
kx = cos(fi*pi/180)*sin(theta*pi/180)*k;
ky = sin(fi*pi/180)*sin(theta*pi/180)*k;
kz = sqrt(k^2 - kx^2 - ky^2);

P_ref_XY = exp(-1i*(kx*Xxy + ky*Yxy));
P_ref_XZ = exp(-1i*(kx*Xxz + kz*Zxz));
P_ref_YZ = exp(-1i*(ky*Yyz + kz*Zyz));

yref = 1.5;
x0 = (-20:dx:20)';
kr = sqrt(k^2-kx^2);

P_xy = zeros(size(Xxy));
P_xz = zeros(size(Xxz));
P_yz = zeros(size(Yyz));
wb = waitbar(0,'Calculating');
for n = 1:length(x0)
   
    waitbar(n/length(x0),wb);
    R_xy = sqrt( ( (Xxy-x0(n)).^2 + Yxy.^2 ) );
    R_xz = sqrt( ( (Xxz-x0(n)).^2 + Zxz.^2 ) );
    R_yz = sqrt( ( x0(n).^2 + Yyz.^2 + Zyz.^2 ) );
    
    P_xy = P_xy + 1/(4*pi)*exp(-1i*kx*x0(n))*exp(-1i*k*R_xy)./R_xy*dx;
    P_xz = P_xz + 1/(4*pi)*exp(-1i*kx*x0(n))*exp(-1i*k*R_xz)./R_xz*dx;
    P_yz = P_yz + 1/(4*pi)*exp(-1i*kx*x0(n))*exp(-1i*k*R_yz)./R_yz*dx;
    
end
close(wb);

Lz = 2*pi/kz;
Ly = 2*pi/ky;
Lr = 2*pi/kr;
A = 1+(Ly/Lz)^2;
B = 2*Ly^2/Lz;
C = Ly^2 - Lr^2;
a = (-B-sqrt(B^2-4*A*C))/(2*A);
%%
pos = [ 0.0 0.05 0.5 1 ;
        0.5 0.05 0.5 0.9 ];
f = figure('Units','points','Position',[200,200,700,330]);

subplot(1,2,1)
mesh(x,y,0*Zxz,real(P_xy),'FaceColor','interp');
set(gca, 'Units','normalized','Position',pos(1,:));
hold on; 
mesh(x,0*Yxy,Zxz,real(P_xz),'FaceColor','interp');
e2 = mesh(0*Yyz,Yyz,Zyz,real(P_yz),'FaceColor','interp','facealpha',.8);
shading interp
axis equal tight
caxis([-1,1]*1e-1)
set(gca,'view',[130,30]);

ol = 15*dx;
xi = [x(1)-ol 0,   0 ]'+dx;
yi = [ 0,  y(1)-ol,  0 ]'+dx;
zi = [ 0,    0,   z(1)-ol  ]'+dx;
Xi = [xi yi zi];
xe = [x(end)+ol, 0,    0 ]';
ye = [ 0,  y(end)+ol,  0 ]';
ze = [ 0,    0,   z(end)+ol  ]';
Xe = [xe ye ze];
dp = Xi-Xe;

arrow3(Xi',Xe','1k')
xlim([xi(1),xe(1)]);
ylim([yi(2),ye(2)]);
zlim([zi(3),ze(3)]);
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');
set(gca,'box','off')
color = get(f,'Color');
set(gca,'XColor',color,'YColor',color,'ZColor',color,'TickDir','out')
grid off

subplot(1,2,2)
pcolor(x,y,real(P_xy));
set(gca, 'Units','normalized','Position',pos(2,:));
axis equal tight
shading interp
caxis([-1,1]*1e-1)
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');
set(gca,'box','off')
color = get(f,'Color');
set(gca,'XColor',color,'YColor',color,'TickDir','out')
grid off
hold on
pcolor_ax =(gca);
cRange = caxis;
P_xy(1:2,:)= 0*P_xy(1:2,:);
[C,~]= contour( x, y , real(P_xy) ,1 , '-k');

headWidth = 5*1.5;
headLength = 5*1.5;
ah = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah,'parent',gca);
set(ah,'position',[x(1)-5*dx, 0, x(end)-x(1)+15*dx, 0]);    
xlim([x(1)-10*dx,x(end)+15*dx])

ah2 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah2,'parent',gca);
set(ah2,'position',[0, y(1)-5*dx, 0, y(end)-y(1)+15*dx]);  
ylim([y(1)-10*dx,y(end)+15*dx])

y0 = 1;
x0 = 0;
r0 = abs(y0);
xstat = r0*kx/kr;
plot(x0,y0,'.k','MarkerSize',20);
plot(-xstat,0,'.k','MarkerSize',20);

ah3 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah3,'parent',gca);
set(ah3,'position',[x0, y0, kx/k, kr/k]);  
ah4 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah4,'parent',gca);
set(ah4,'position',[x0, y0, 0, kr/k]);  
ah5 = annotation('arrow',...
    'headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth);
set(ah5,'parent',gca);
set(ah5,'position',[x0, y0, kx/k, 0]);  


line([-xstat,0],[x0,y0],'Color','black','LineStyle','--')

ylim([y(1)-10*dx,y(end)+15*dx])

set(gcf,'PaperPositionMode','auto');
print( '-r300', fullfile( '../..','Figures/High_freq_approximations','greens_stat_pos_2' ) ,'-dpng')