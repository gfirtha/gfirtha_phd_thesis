function [ p,xp,yp ] = make_path( x_a, y_a, span )

D = sqrt(diff(x_a).^2 + diff(y_a).^2);

xp = 0;
yp = 0;
for n = 1 : length(x_a)-1
    xp = [ xp(1:end-1) linspace( x_a(n),x_a(n+1), round(D(n)*100) ) ];
    yp = [ yp(1:end-1) linspace( y_a(n),y_a(n+1), round(D(n)*100) ) ];
end

xp = smooth(xp,span);
yp = smooth(yp,span);
p = [0;cumsum( sqrt(diff(xp).^2 + diff(yp).^2) )];

end

