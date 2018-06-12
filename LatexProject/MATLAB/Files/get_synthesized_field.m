function [ field_synth ] = get_synthesized_field( t,s, Amp, xs, antialiasing_prefilter ,loudspeaker_type, reproduction_filter, Dref,kn,r0, xrec0, x_ssd, N_ssd, X,Y, kv, kv0 )

c = 343.1;
fs =  1/mean(diff(t));
ds = max(sqrt(sum(diff(x_ssd,1).^2,2)));

Nt = length(t);
w = fftshift( 2*pi*(-Nt/2:Nt/2-1)'/(Nt)*fs );

H = sqrt(8*pi*1i*w/c);
h = fftshift(ifft(H));

s_filt = conv(s,h,'same');
Dx0 = Amp/(4*pi)*sqrt(Dref).*kn./r0;

%%
fi = (-pi:pi/180:pi);
[FI,W] = meshgrid(fi,w);
% Baffled circular piston model
if strcmp(loudspeaker_type,'piston')
    R_ls = ds/2;
    kRsinFi = bsxfun( @times, w/c*R_ls , sin(fi));
    Dir_char_mx = 2*besselj(1,kRsinFi)./(kRsinFi);
    Dir_char_mx(isnan(Dir_char_mx)) = 1;
    
    % Ideal reproduction loudspeakers
elseif strcmp(loudspeaker_type,'ideal')
    F1 = zeros(length(w),length(fi));
    KX = (W/c).*sin(FI);
    Dir_char_mx = get_transfer(KX,pi/ds,0.5);
end
%
if strcmp(reproduction_filter,'on')
dir_ir = ifftshift( real( ifft(  Dir_char_mx,[],1 ) ) );
Limp = 500;
dir_ir = dir_ir( round(end/2)-(Limp/2-1):round(end/2)+Limp/2,:).*tukeywin(Limp,0.5);
s_filt_ls = zeros(length(t),length(fi));
for n = 1 : length(fi)
    s_filt_ls(:,n) = conv(s_filt,dir_ir(:,n),'same');
end
elseif strcmp(reproduction_filter,'off')
   s_filt_ls = repmat(s_filt,1,length(fi)); 
end

if strcmp(antialiasing_prefilter,'on')
% Anti aliasing filtering
w_c = (pi/ds)*c./(abs(kv-kv0));
w_c(w_c>2*pi*fs/2 *0.99) = 2*pi*fs/2 *0.99;
s_filt_ls_aa = zeros( size(s_filt_ls,1),size(s_filt_ls,2), length(Dx0) );

Nfil = 20;
filtz = zeros(Nfil+1,length(x_ssd));
for n = 1 : length(x_ssd)
    f = design_filter(Nfil,0.75*w_c(n)/(2*pi),1*w_c(n)/(2*pi),fs);
    filtz(:,n) = f.Numerator';
    for m = 1 : size(s_filt_ls,2)
        s_filt_ls_aa (:,m,n) = conv(s_filt_ls(:,m), filtz(:,n) ,'same' );
    end
end
elseif strcmp(antialiasing_prefilter,'off')
    s_filt_ls_aa = repmat(s_filt_ls,1,1,length(x_ssd));
end
%%
t_0 = norm(xrec0-xs,2)/c;

field_synth = zeros(size(X));
tic
[FI,T] = meshgrid(fi,t);
wb = waitbar(0,'Calculating');
for n = 1 : length(Dx0)
    
    waitbar(n/length(Dx0),wb);
    S = squeeze( s_filt_ls_aa(:,:,n) );
    R = sqrt( (X-x_ssd(n,1)).^2 + (Y-x_ssd(n,2)).^2 );
    
    xr = bsxfun( @times, [X(:)-x_ssd(n,1), Y(:)-x_ssd(n,2)] , 1./sqrt( (X(:)-x_ssd(n,1)).^2 + ( Y(:)-x_ssd(n,2) ) .^2 ) );
    fi_field = reshape( acos( xr*N_ssd(n,:)' ), size(X,1), size(X,2) );

    [FI2,T] = meshgrid(fi,t');
    t_ret = t_0 - R/c - r0(n)/c;
    S0  =  interp2(FI2,T,S,fi_field,t_ret,'linear');
    field_synth = field_synth + 1/(4*pi)*Dx0(n).*S0 ./R*ds;
end
close(wb)
toc



end

