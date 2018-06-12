function [ d_wfs_aa ] = anti_aliasing_fd(fs,x0,d_wfs,ts,wc,wlen,hop,nfft)

d_wfs_aa = zeros(size(d_wfs));
wb = waitbar(0,'Antialiasing filtering in the STFT domain');
for n = 1 : length(x0)
    waitbar(n/length(x0),wb);
    [Dxw,f_stft,t_stft] = stft(d_wfs(:,n), wlen, hop, nfft,fs);
    fc = interp1(ts,wc(:,n),t_stft,'linear','extrap')/(2*pi);
    [~,w_ind ]= min(abs(bsxfun(@minus, f_stft', fc)),[],1);
    AAF_mx = zeros(size(Dxw));
    for i = 1 : size(AAF_mx,2)
        temp = fftshift(tukeywin(w_ind(i)*2,1));
        AAF_mx(1:w_ind(i),i) = temp(1:end/2);
    end
    Dfxw = Dxw.*AAF_mx;
    Dfilt = istft(Dfxw,wlen, hop, nfft,fs);
    d_wfs_aa(1:length(Dfilt),n) = Dfilt;
end
close(wb);


end

