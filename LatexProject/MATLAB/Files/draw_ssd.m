function [  ] = draw_ssd( f, x0, n0, size )

for n = 1 : length(x0)
draw_ls(f,[x0(n,1),x0(n,2)], [n0(n,1),n0(n,2)]', size)
end


end

