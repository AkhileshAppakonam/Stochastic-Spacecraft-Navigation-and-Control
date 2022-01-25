function [dU] = gradU(x,y,z,C20,C22,mu)
% U = mu/r - (mu*C20*(x^2+y^2 - 2*z^2))/(2*r^5) + (3*mu*C22*(x^2-y^2)/(r^5));

rn = x^2 + y^2 + z^2;
Ux = (6*C22*mu*x)/rn^(5/2) - (C20*mu*x)/rn^(5/2) - (mu*x)/rn^(3/2) + ...
     (5*C20*mu*x*(x^2 + y^2 - 2*z^2))/(2*rn^(7/2)) - (15*C22*mu*x*(x^2 - y^2))/rn^(7/2);
 
Uy = (5*C20*mu*y*(x^2 + y^2 - 2*z^2))/(2*rn^(7/2)) - (C20*mu*y)/rn^(5/2) - ...
     (6*C22*mu*y)/rn^(5/2) - (mu*y)/rn^(3/2) - (15*C22*mu*y*(x^2 - y^2))/rn^(7/2);
 
Uz = (2*C20*mu*z)/rn^(5/2) - (mu*z)/rn^(3/2) + ...
     (5*C20*mu*z*(x^2 + y^2 - 2*z^2))/(2*rn^(7/2)) - (15*C22*mu*z*(x^2 - y^2))/rn^(7/2);

% dU = [Ux;Uy;Uz];


% POINT MASS
% U = mu/r 
dU = [-(mu*x)/rn^3; -(mu*y)/rn^3; -(mu*z)/rn^3];
end

