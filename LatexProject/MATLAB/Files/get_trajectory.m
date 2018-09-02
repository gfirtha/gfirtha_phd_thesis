function [ Xs ] = get_trajectory( p,xp,yp, t, v )

t0    = zeros(length(t),1);
t0(1) = 0;

h = mean(diff(t));
if numel(v) ~= 1
    wb = waitbar(0,'Calculating  source location vector');
    for n = 1 : length(t)-1
        waitbar(n/length(t),wb);
        t0(n+1) = t0(n) + interp1(t,v,t0(n),'linear','extrap')*h;
    end
    close(wb)
    xs = interp1(p,xp,t0,'linear','extrap');
    ys = interp1(p,yp,t0,'linear','extrap');
else
    
    xs = interp1(p,xp,t*v,'linear','extrap');
    ys = interp1(p,yp,t*v,'linear','extrap');
end
Xs = [xs ys];


end

