function [  ] = draw_ls( f, center, normal, size )

normal = normal/norm(normal);

fi = -atan2d(normal(1)*1,normal(2)*1)*pi/180;
x1 = [   -1    1   1   -1 ]'*size ;
y1 = [  -1.8  -1.8  -1  -1 ]'*size ;
x1_ = cos(fi)*x1 - sin(fi)*y1;
y1_ = sin(fi)*x1 + cos(fi)*y1;

x2 = [   -1    1  1.5 -1.5 ]'*size ;
y2 = [   -1 -1    0   0 ]'*size ;
x2_ = cos(fi)*x2 - sin(fi)*y2;
y2_ = sin(fi)*x2 + cos(fi)*y2;

%figure(f);
hold on
fill(  x1_ + center(1),y1_ + center(2),ones(1,3)*0.2 )
fill(  x2_ + center(1),y2_ + center(2),ones(1,3)*0.5 )
end

