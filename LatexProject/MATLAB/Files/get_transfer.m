function [ H ] = get_transfer( in, lim, L )
H = cos((in./lim-L)/(1-L)*pi/2).^2;
H( abs(in./lim)< L ) = 1;
H( abs(in./lim)>= 1 ) = 0;


end

