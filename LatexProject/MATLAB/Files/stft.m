function [ out, t_ ] = stft( input, w, fs, W, step )

nW = round(W*fs);
stepN = round(step*fs);
N = floor( (length(input) - nW + stepN)/stepN ) ;

t = (0:nW-1)'/fs;
tr_mx = exp( - 1i* w * t' )*1/fs;
out = zeros( size(w,1), N );

for n = 1:N 
    in = input( (n-1)*stepN + 1 : (n-1)*stepN + nW ).';
    out(:,n) = tr_mx*(in.*tukeywin(nW,1));
end
t_ = W/2+(0:N-1)*step;

end

