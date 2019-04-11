clear
close all


fi = linspace(0,30,30);
fi0 = 30;

C = ( tand(fi)/tand(fi0) );
A1 = ones(size(fi));
A2 = A1.*(1-C)./(C+1);
%
Avec = [A1;A2];
Avec = bsxfun(@times, Avec, 1./sqrt(sum(Avec.^2,1)));
d = 20*log10(Avec(1,:)./Avec(2,:));

d0 = linspace(0,20,30);

D_leakey_I = [0,3.8,8,12.5,17]';
Leakey_I = [0,7.5,15,22,28.5]';
D_Bl_I = [0,3.5,7.75,11.25,14]';
BrittainLeakey_I = [0,10,20,28,30]';
D_deboer_I = [0,4.5,11,14];
DeBoer_I = [0,10,23,27]';

T_deboer = [0 10   20  30   40    50  60 0.2*343]/343*10;
DeBoer_T = [0, 4, 7.25, 10.5, 13.5, 16, 18.5,  20];

T_Leakey = [ 0, 0.3, 0.65, 1.2, 1.75 ]*0.5;
Leakey_T = [0 7.5, 15, 23,27];
%%
ftsize = 9;

f = figure('Units','points','Position',[200,200,204,214]);

set(f,'defaulttextinterpreter','latex')

pos = [ 0.18 0.63 0.75 .35
        0.18 0.11 0.75 .35];
p1 = axes('Units','normalized','Position',pos(1,:));
plot(d,fi,'k','LineWidth',1)
hold on
plot(linspace(0,D_Bl_I(end),30), spline(D_Bl_I, BrittainLeakey_I,linspace(0,D_Bl_I(end),30)),':k','LineWidth',1);
plot(linspace(0,D_leakey_I(end),30), spline(D_leakey_I, Leakey_I,linspace(0,D_leakey_I(end),30)),'--k','LineWidth',1);
plot(linspace(0,D_deboer_I(end),30), spline(D_deboer_I, DeBoer_I,linspace(0,D_deboer_I(end),30)),'-.k','LineWidth',1);
xlim([0,20])    
ylim([0,30])
grid on
xlabel( 'Interchannel amplitude difference $[\mathrm{dB}]$' , 'FontSize', ftsize );
ylabel( '\parbox{.65in}{Phantom direction}$\phi_p$ [$^{\circ}$]' , 'FontSize', ftsize );
legend('Tangent law','Brittain','Leakey','de Boer','Location','SouthEast')
set(gca,'FontName','Times New Roman');

p2 = axes('Units','normalized','Position',pos(2,:));
plot(linspace(0,T_Leakey(end),30), spline(T_Leakey, Leakey_T,linspace(0,T_Leakey(end),30)),'--k','LineWidth',1);
hold on
plot(linspace(0,T_deboer(end),30), spline(T_deboer, DeBoer_T,linspace(0,T_deboer(end),30)),'-.k','LineWidth',1);
%xlim([0,20])    
ylim([0,30])
grid on
xlabel( 'Interchannel time difference$ [\mathrm{ms}]$' , 'FontSize', ftsize );
ylabel( '\parbox{.65in}{Phantom direction}$\phi_p$ [$^{\circ}$]' , 'FontSize', ftsize );
legend('Leakey','de Boer','Location','SouthEast')
set(gca,'FontName','Times New Roman');


allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300',fullfile( '../..','Figures/Introduction','stereo_b' ) ,'-dpng')

