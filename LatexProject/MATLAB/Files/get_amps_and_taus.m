function [ A, Tau, wc ] = get_amps_and_taus( ts, x0,n0,v0, xs, Tau0,c,Rref )
% The function calculates the amplitude factors and delay calues for the
% moving source driving functions at time instants ts


Tau = zeros(length(ts),length(x0));
A     = zeros(length(ts),length(x0));
wc    = zeros(length(ts),length(x0));

Tau_n = Tau0;
dx0 = mean(sqrt(sum(diff(x0,1).^2,2)));
fs0 = 1/mean(diff(ts));

wb = waitbar(0,'Calculating amplitude and propagation delays');
for n = 1 : length(ts)
    waitbar(n/length(ts),wb);
    
    xs_t = interp1( ts,xs(:,1), ts(n)-Tau_n, 'linear','extrap' );
    ys_t = interp1( ts,xs(:,2), ts(n)-Tau_n, 'linear','extrap' );
    vx = (xs_t - interp1( ts,xs(:,1), ts(n)-1/fs0-Tau_n, 'linear','extrap' ))*fs0;
    vy = (ys_t - interp1( ts,xs(:,2), ts(n)-1/fs0-Tau_n, 'linear','extrap' ))*fs0;
    Dvx = x0-[xs_t ys_t];
    
    R = sqrt( sum( Dvx.^2,2) );
    Vv = 1/c*sum([vx vy].*(Dvx),2);
    Delta = R - Vv;
    
    D = xs_t.*x0(:,2) - ys_t.*x0(:,1);
    
    xref = ( D.*Dvx(:,2) - sign(Dvx(:,2)).*Dvx(:,1).*sqrt(Rref^2.*R.^2-D.^2 ))./R.^2;
    yref = (-D.*Dvx(:,1) - abs(Dvx(:,2)).*sqrt(Rref^2.*R.^2-D.^2 ))./R.^2;
    dref = sqrt( sum( (x0-[xref yref]).^2, 2)  );
    %   dref = Rssd;
    
    Kn = sum(n0.*Dvx ,2);
    Kn =  Kn.*(Kn>0);
    A(n,:) = sqrt(dref./(dref+R)).* Kn./Delta.^(3/2)*dx0;
    
    Tau(n,:) = Tau_n;
    Tau_n = Tau_n - 1/fs0*Vv./Delta;
    
    kv = sum(bsxfun(@times, Dvx, 1./R).*v0,2);
    wc(n,:) = pi/dx0*c./abs(kv);
    
end
close(wb);

end

