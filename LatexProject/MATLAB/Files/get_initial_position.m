function [ Tau0 ] = get_initial_position( v,c, x_a, x0 )

M = v/c;
phi0 = acos( (x_a(1,2)-x_a(1,1))./sqrt( (x_a(1,2)-x_a(1,1)).^2 + (x_a(2,2)-x_a(2,1)).^2  ) );
x_ =  cos(phi0)*(x0(:,1)-x_a(1,1)) + sin(phi0)*(x0(:,2)-x_a(2,1));
y_ = -sin(phi0)*(x0(:,1)-x_a(1,1)) + cos(phi0)*(x0(:,2)-x_a(2,1));
Tau0 = (M*x_+sqrt( x_.^2 + y_.^2*(1-M^2)))./(c*(1-M^2));


end

