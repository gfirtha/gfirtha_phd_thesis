function [ x_ssd, N_ssd, ds ] = create_equidist_SSD( ssd_contour, No_ls )

x_sec_contour = ssd_contour(:,1);
y_sec_contour = ssd_contour(:,2);
Dx = x_sec_contour(2:end)-x_sec_contour(1:end-1);
Dy = y_sec_contour(2:end)-y_sec_contour(1:end-1);
ArcLength = sqrt( Dx.^2 + Dy.^2);
CurveLength = cumsum ( ArcLength );

rep_elem = find(( Dx == 0).*( Dy == 0));
x_ssd_temp = x_sec_contour;
x_ssd_temp( rep_elem ) = [];
y_ssd_temp = y_sec_contour;
y_ssd_temp( rep_elem )  = [];
CurveLength( rep_elem )  = [];

x_new_nodes = interp1( [0; CurveLength], x_ssd_temp, linspace(0,CurveLength(end),No_ls+1) )';
y_new_nodes = interp1( [0; CurveLength], y_ssd_temp, linspace(0,CurveLength(end),No_ls+1) )';

x_ssd = [(x_new_nodes(1:end-1)+x_new_nodes(2:end))/2 (y_new_nodes(1:end-1)+y_new_nodes(2:end))/2];

DX = diff(x_new_nodes);   DY = diff(y_new_nodes);
%
d_ssd = bsxfun(@times,  [DX DY] , 1./sqrt(DX.^2+DY.^2));
N_ssd = ( [0 -1; 1 0] *d_ssd')';

ds = CurveLength(end)/No_ls;

end