function [a,e,i_d,Omega_d,omega_d,M_d,E_d,theta_d,l_d,Pi_d,u_d]= cart2kepl(x,y,z,xdot,ydot,zdot,mu)
% Inputs
% X=[x,y,z]                       position vector [m]
% V=[xdot,ydot,zdot]              velocity vector [m/s]
% mu                              planet gravitational parameter [m^3 s^-2]

% Outputs
% a                               semi?major axis [m]
% e                               eccentricity [?]
% i_d                             inclination [deg]
% Omega_d                         longitude of the ascending node (RAAN) [deg]
% omega_d                         argument of perigee [deg]
% M_d                             meananomaly[deg]
% l_d                             true longitude [deg] 
% Pi_d                            longitude of perigee [deg]
% u_d                             argument of latitude [deg]


format long g

% introduce unit vectors
i_ =[1;0;0]; % unit vector in the X direction
j_ =[0;1;0]; % unit vector in the Y direction
k_ =[0;0;1]; % unit vector in the Z direction

% introduce position vector [m] , velocity vector [m/s ]
r_ =[x;y;z]; % position vector
V_ =[xdot;ydot;zdot];  % velocity vector
r=norm(r_) ;  % norm of position vector
V=norm(V_);  % norm of velocity vector

% compute angular momentum vector [m?2/s ] 
H_=cross(r_,V_); % angular momentum vector
H=norm(H_) ; % angular momentum norm

% compute semi?major axis [m] 
a=r/(2-r*V^2/mu) ;

% compute eccentricity vector
e_ =1/mu*((V^2-mu/r)*r_-dot(r_,V_)*V_); % eccentricity vector
e=norm(e_); % eccentricity norm

% compute inclination
i=acos(dot(k_ ,H_ )/H);

% compute ascending node vector
n_ =cross(k_ ,H_ ); % ascending node vector 
n=norm(n_ ); % ascending node norm

% compute RAAN [rad]
Omega=acos(dot(i_ ,n_ )/n); 
if dot(n_ ,j_ )<0
Omega=2* pi -Omega ; % quadrant check
end

% compute argument of perigee [ rad ]
omega=acos(dot(n_ ,e_ )/(n*e));
if dot(e_ ,k_ )<0    %quadrant check
omega=2*pi-omega ;
end

% compute true anomaly [ rad ]
theta=acos(dot(e_ ,r_ )/(e*r));
if dot(r_ ,V_ )<0
    theta=2*pi-theta ; %quadrant check
end

% compute eccentric anomaly [ rad ]
E= asin( 1/e * sqrt( 1/( mu * a ) ) * ( x * xdot + y * ydot + z * zdot ) ) ; 
if 1/e*sqrt(1/(mu*a))*(x*xdot+y*ydot+z*zdot)>0
   if (1-r/a)/e<0
       E=pi-E;
   end
elseif 1/e*sqrt (1/(mu*a))*(x*xdot+y*ydot+z*zdot)<0 % quadrant check
if (1-r/a)/e<0 
    E=pi-E;
else (1-r/a)/e>0; 
    E=2*pi+E;
end
end

% compute mean anomaly [ rad ] 
M=E-e*sin(E); % mean anomaly relation

% compute true longitude [ rad ]
l = acos(r_(1)/r);
idx = r_(2) < 0; 
if any(idx); 
    l(idx) = 2*pi - l(idx); 
end
% l=acos(dot(i_ ,r_ )/r); 
% if dot(i_ ,r_ )<0
% l=2*pi-l ; % quadrant check end
% end

% compute longitute of perigee [ rad ]
Pi=acos(dot(e_ ,i_ )/e); 
if dot(e_ ,i_ )<0
Pi=2*pi-Pi ; % quadrant check end
end

% compute argument of latitude [ rad ]
u=acos(dot(n_ ,r_ )/(n*r)); 
if dot(n_ ,r_ )<0
u=2*pi-u; % quadrant check end
end


% Account for cases when orbit is near?circular and/or near?equatorial 
tol=1E-8; % define a limit for considering a orbit near?circular and/or near?equatorial
if e<tol && mod(i ,pi)<tol % orbit is near?circular and near?equatorial
Omega=NaN ; 
omega =NaN ; 
theta=l ;
end
if e<tol && mod(i , pi)>tol % orbit is near?circular 
omega =NaN ;
theta=u;
end
if e>tol && mod(i ,pi)<tol % orbit is near?equatorial 
Omega=NaN ;
omega=Pi ;
end
% Convert computed angle from radians to degree
E_d=rad2deg(E);
i_d=rad2deg(i);
theta_d=rad2deg(theta);
M_d=rad2deg (M) ;
Omega_d=rad2deg (Omega) ; 
omega_d=rad2deg(omega);
l_d=rad2deg(l); 
Pi_d=rad2deg(Pi); 
u_d=rad2deg(u); 
end



