function [F,M] = Pan_GravityModel(r,R,V,t)

global II mu Wp C22 C20

PRA = wrapTo2Pi(norm(Wp) *t);
%[rad]Earth rotation angle at mission commencement.

PcpfToPci = [cos(PRA), -sin(PRA), 0; sin(PRA), cos(PRA), 0; 0, 0, 1];
%[]Matrix needed to transform the current position from E.C.E.F. coordinates to E.C.I. coordinates.

% Data
J = II(1:3,1:3);
m = II(4,4);
rn = norm(r);
JJ = 1/2*trace(J)*eye(3) + J;
R_bar = R*PcpfToPci;
pc = R'*r;
rc = PcpfToPci'*r;
x = rc(1);
y = rc(2);
z = rc(3);

% a_sat = + 2*cross(Wp,V(4:6) - cross(Wp,pc)) + cross(Wp,cross(Wp,pc));
%[N]Tertiary gravity force due to the Coriolis force and Saturn.

F2gac =  m * [ -mu*C20*x/rn^5 + 5*mu*C20*x*(x^2+y^2-2*z^2)/(2*rn^7) + 6*mu*C22*x/rn^5 - 15*mu*C22*x*(x^2-y^2)/rn^7;
               -mu*C20*y/rn^5 + 5*mu*C20*y*(x^2+y^2-2*z^2)/(2*rn^7) - 6*mu*C22*y/rn^5 - 15*mu*C22*y*(x^2-y^2)/rn^7;
              2*mu*C20*z/rn^5 + 5*mu*C20*z*(x^2+y^2-2*z^2)/(2*rn^7) - 15*mu*C22*z*(x^2-y^2)/rn^7 - mu*z/rn^3]; %- mu*z/rn^3
%[N]Secondary gravity force due to spherical harmonics
          
Fpgac = -m*mu/rn^3 *pc -3*(mu/rn^5)*JJ*pc + 15/2*(mu/rn^7 * pc' * J * pc)*pc;
%[N]Primary gravity force

F = Fpgac + R_bar*F2gac;%- m*a_sat;
%[N]Sum of all forces rotated into Pan Centered Pan Fixed Frame

M = 3 *(mu/rn^5)* cross(pc,J*pc);
%[]Current Moment of Inertial

end

