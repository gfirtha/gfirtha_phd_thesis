function [ Dx_filt ] = anti_aliasing_td( Dx, t, ts,wc0, Nfil, step )

fs = 1/mean(diff(t));
for n = 1 : size(wc0,2)
    wc(:,n) = interp1(ts,wc0(:,n),t,'linear','extrap');
end

fc = wc/(2*pi);
fc(fc>0.95*fs/2) = 0.95*fs/2;

zi = [];
out = zeros(size(Dx,1),1);
N = floor(size(Dx,1)/step);
Dx_filt = zeros(size(Dx));
wb = waitbar(0,'Antialiasing filtering using FIR filters');
for i = 1 : size(Dx,2)
    waitbar(i/size(Dx,2),wb);
    for n = 1:N
        Fc = mean(fc((n-1)*step+1:n*step,i));
     %   fir = fir1(Nfil,Fc);
        fir = design_filter(Nfil,0.75*Fc,1*Fc,fs);
        in = Dx( (n-1)*step+1:n*step,i )';
        [out( (n-1)*step+1:n*step ),zi] = filter(fir.Numerator',1,in,zi);
    end
    Dx_filt(:,i) =  circshift(out,0);
end
close(wb)


end

