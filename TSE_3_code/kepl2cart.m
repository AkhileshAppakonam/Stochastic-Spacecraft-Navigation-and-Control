function [x,y,z,xdot,ydot,zdot] = kepl2cart(a,e,i_d,Omega_d,omega_d,M_d,l_d,Pi_d,u_d,mu)

% Inputs
% a                               semi?major axis [m]
% e                               eccentricity [-]
% i_d                             inclination [deg]
% Omega_d                         longitude of the ascending node (RAAN) [deg]
% omega_d                         argument of perigee [deg]
% M_d                             meananomaly[deg]
% l_d                             true longitude [deg] 
% Pi_d                            longitude of perigee [deg]
% u_d                             argument of latitude [deg]
% mu                              planet gravitational parameter [m^3 s^-2]

% Outputs
% X=[x,y,z]                       position vector [m]
% V=[xdot,ydot,zdot]              velocity vector [m/s]

H= sqrt( a * ( 1 - e^2 ) * mu ) ;
% convert angles to radians
i=deg2rad(i_d);
Omega=deg2rad(Omega_d) ; 
omega=deg2rad(omega_d);
M=deg2rad(M_d); 
l=deg2rad(l_d);
Pi=deg2rad(Pi_d);
u=deg2rad(u_d);

% compute eccentric anomaly E and true anomaly theta from mean anomaly M
% syms E theta
E0 = 0;
Tol=1;
while  Tol> 1e-6
    E1 = M + e*sin(E0);
    Tol = abs(E1-E0);
    E0=E1;
end
E = E1;
theta= 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); %relation between E and theta

% E_eq=M==E-e*sin(E); % mean anomaly relation
% theta_eq=tan(theta/2)==sqrt((1+e)/(1-e))*tan(E/2); %relation between E and theta
% sol=vpasolve([E_eq, theta_eq],[E,theta]);
% E= double( sol.E ) ;
% theta = double(sol.theta);

% Account for cases when orbit is near?circular and/or near?equatorial 
tol=1E-8; % define a limit for considering a orbit near?circular and/or near?equatorial
if e<tol && mod(i ,pi)<tol % orbit is near?circular and near?equatorial
Omega=0; omega=0; theta=l ;
end
if e<tol && mod(i , pi)>tol % orbit is near?circular omega=0;
theta=u;
end
if e>tol && mod(i ,pi)<tol % orbit is near?equatorial Omega=0;
omega=Pi ;
end

% compute radius vector r
r=(a*(1-e^2))/(1+e*cos(theta)) ;

% compute the coordinates in the rectangular non?rotating reference frame xi , eta , zeta
xi=r*cos(theta); eta=r*sin(theta);

% compute the velocities in the rectangular non?rotating reference frame xi , eta , zeta
xi_dot=-mu/H*sin(theta); eta_dot=mu/H*(e+cos(theta)); 

% compute the parameters for the rotation matrix
l1 = cos ( omega ) * cos ( Omega )- sin ( omega ) * sin ( Omega ) * cos ( i );
m1= cos ( omega ) * sin ( Omega ) + sin ( omega ) * cos ( Omega ) * cos ( i ) ; 
n1=sin (omega)*sin ( i ) ;
l2 =- sin ( omega ) * cos ( Omega ) - cos ( omega ) * sin ( Omega ) * cos ( i ) ; 
m2=- sin ( omega ) * sin ( Omega ) + cos ( omega ) * cos ( Omega ) * cos ( i ) ; 
n2= cos ( omega ) * sin ( i ) ;

% build the rotation matrix RR
RR=[l1,l2; m1,m2; n1,n2];

% do the actual conversion
X=RR*[xi;eta]; 
V=RR*[xi_dot;eta_dot];

% define the output components
x=X(1) ; y=X(2) ; z=X(3) ; xdot=V(1) ; ydot=V(2) ; zdot=V(3) ;
end

