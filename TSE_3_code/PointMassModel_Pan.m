%===================================================================================================
%[]FUNCTION NAME: NSER2BPModel_Kalman.m
%[]AUTHOR: Matthew Wittal
%[]CREATED: 12/09/2019
%[]REVISED: --/--/----
%===================================================================================================
%[]FUNCTION DESCRIPTION:
%This function evaluates the equations of motion for the relative two-body problem while taking
%Earth rotation and spherical harmonics into account.
%===================================================================================================
%[]INPUT VARIABLE:
%(t)|Current time vector.
%---------------------------------------------------------------------------------------------------
%(S)|Current state vector.
%===================================================================================================
%[]OUTPUT VARIABLES:
%(dS)|Current differential state vector
%===================================================================================================
%[]VARIABLE FORMAT:
%(t)|Scalar {1 x 1}.
%---------------------------------------------------------------------------------------------------
%(S)|Column Vector {6 x 1}.
%---------------------------------------------------------------------------------------------------
%(dS)|Column Vector {6 x 1}
%===================================================================================================
%[]AUXILIARY FUNCTIONS:
%None.
%===================================================================================================
%[]COMMENTS:
%None.
%===================================================================================================
function dS = PointMassModel_Pan(t,S)
global II m mu Wp C20 C22

PRA = wrapTo2Pi(norm(Wp) *t);
%[rad]Earth rotation angle at mission commencement.

PcpfToPci = [cos(PRA), -sin(PRA), 0; sin(PRA), cos(PRA), 0; 0, 0, 1];
%[]Matrix needed to transform the current position from E.C.E.F. coordinates to E.C.I. coordinates.

r = S(1:3);

V = [0;0;0; S(4:5)];

R = eye(3);

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
 
% a_sat = + 2*cross(Wp,S(4:6) - cross(Wp,pc)) + cross(Wp,cross(Wp,pc));
%[N]Tertiary gravity force due to the Coriolis force and Saturn.

F2gac =  m * [ -mu*C20*x/rn^5 + 5*mu*C20*x*(x^2+y^2-2*z^2)/(2*rn^7) + 6*mu*C22*x/rn^5 - 15*mu*C22*x*(x^2-y^2)/rn^7;
               -mu*C20*y/rn^5 + 5*mu*C20*y*(x^2+y^2-2*z^2)/(2*rn^7) - 6*mu*C22*y/rn^5 - 15*mu*C22*y*(x^2-y^2)/rn^7;
              2*mu*C20*z/rn^5 + 5*mu*C20*z*(x^2+y^2-2*z^2)/(2*rn^7) - 15*mu*C22*z*(x^2-y^2)/rn^7 - mu*z/rn^3 ]; %- mu*z/rn^3
%[N]Secondary gravity force due to spherical harmonics
          
Fpgac = -m*mu/rn^3 *pc -3*(mu/rn^5)*JJ*pc + 15/2*(mu/rn^7 * pc' * J * pc)*pc;
%[N]Primary gravity force

F =  Fpgac + R_bar*F2gac;
%[N]Sum of all forces rotated into Pan Centered Pan Fixed Frame
    
dS(1:3)= S(4:6);
% dS(4:6)= -mu/rn^3 * r;
dS(4:6)= F/m;%- a_sat;
dS=dS';
end
%===================================================================================================