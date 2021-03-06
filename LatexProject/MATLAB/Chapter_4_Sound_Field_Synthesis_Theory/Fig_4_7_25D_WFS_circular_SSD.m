clear
close all

f = 1e3;
speed_of_sound = 343.1;
omega = 2*pi*f;
k = omega/speed_of_sound;

R0 = 2;
dfi = 2*pi/1000;
fi = (0:dfi:2*pi-dfi)';
x0 = R0*cos(fi);
y0 = R0*sin(fi);
xssd = [x0 y0];
n_ssd = [ -cos(fi) -sin(fi) ];

x = linspace(-3,3,200);
y = linspace(-3,3,200);
[X,Y] = meshgrid(x,y);

% Reference a plane wave to a circle
alpha_pw = 0*pi/180;  
K = [cos(alpha_pw); sin(alpha_pw)]*k;
kn = n_ssd * K;
a0 = double((kn)>0);
Kh = K/k;

Rref = 1.5;
dref = R0* ( kn/k - sqrt( (kn/k).^2 + (Rref/R0)^2 - 1 ) );

Dx0 = a0.*sqrt(8*pi/(1i*k)).*sqrt(dref).*1i.*kn.*exp(-1i*( K(1)*x0 + K(2)*y0 ));
%%
field_pw = zeros(size(X));
field_ps = zeros(size(X));

wb = waitbar(0,'Calculating radiated field');
for n = 1:length(x0)
   
    waitbar(n/length(x0),wb);
    r_ssd = sqrt(( X - x0(n) ).^2 + ( Y - y0(n) ).^2 );
    field_pw = field_pw + R0/(4*pi)*Dx0(n)*exp(-1i*k*r_ssd)./r_ssd*dfi;                              
    
end
close(wb);

field_pw_ref = exp(-1i*k*(cos(alpha_pw)*X+sin(alpha_pw)*Y));

%%
%x_pw = x0 + dref*cos(alpha_pw);
%y_pw = y0 + dref*sin(alpha_pw);

x_pw = Rref*cos(fi((kn)>=0));
y_pw = Rref*sin(fi((kn)>=0));
x_pw2 = Rref*cos(fi((kn)<0));
y_pw2 = Rref*sin(fi((kn)<0));

[~,i] = max(sqrt(sum(diff(x_pw2).^2,2)));
x_pw2 = [x_pw2(i+1:end,:);x_pw2(1:i,:)];
y_pw2 = [y_pw2(i+1:end,:);y_pw2(1:i,:)];

ftsize = 9;
f = figure('Units','points','Position',[200,120,407,188]);
set(f,'defaulttextinterpreter','latex')

pos = [ 0.045 0.16 0.4 .76 
        0.52  0.16 .5 .76  ];


p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(field_pw));
shading interp
axis equal tight
hold on
plot(x0,y0,'k','LineWidth',1)
plot(x_pw,y_pw,'--w','LineWidth',1)
plot(x_pw2,y_pw2,':w','LineWidth',1)
caxis([-1,1])
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);

p2 = axes('Units','normalized','Position',pos(2,:));
pcolor(x,y,20*log10(abs(field_pw_ref-field_pw)));
shading interp
axis equal tight
hold on
plot(x0,y0,'k','LineWidth',1)
plot(x_pw,y_pw,'--w','LineWidth',1)
plot(x_pw2,y_pw2,':w','LineWidth',1)
caxis([-25,20])
xlabel( '$x$ [m]', 'FontSize', ftsize );
ylabel( '$y$ [m]', 'FontSize', ftsize );
%title('$20\mathrm{log}_{10}\left( P_{\mathrm{synth,pw}}(\mathbf{x},\omega)-P_{\mathrm{ref,pw}}(\mathbf{x},\omega) \right)$'...
%                                                , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);

xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
c = colorbar;
title(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize);

set(gcf,'PaperPositionMode','auto');
print(  '-r300', fullfile( '../..','Figures/SFS_theory','25D_WFS_circular_SSD' ) ,'-dpng')