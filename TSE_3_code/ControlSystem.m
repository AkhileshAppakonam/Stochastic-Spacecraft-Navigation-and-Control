function [psi,u] = ControlSystem(g,V)

% -------- Variables --------
global target_orbit e kappa a II K1 K2 controller target_attitude horizon
R = g(1:3,1:3);
r = g(1:3,4);
u = zeros(6,1);




% -------- Guidance System --------
d = vecnorm(target_orbit(:,1:3)'-r);
[M,I] = min(d);

if M<5
k=horizon;
else
k=0;
end

if I+horizon > length(target_orbit)  
% TrajRef=target_orbit(1,1:3)';
% VelRef=target_orbit(1,4:6)';
target_orbit = [target_orbit;target_orbit];
TrajRef=target_orbit(I+k,1:3)';
VelRef=target_orbit(I+k,4:6)';
else
TrajRef=target_orbit(I+k,1:3)';
VelRef=target_orbit(I+k,4:6)';
end

% TrajRef=[-125;70;-20];
% VelRef= zeros(3,1);

% figure(101)
% plot3(TrajRef(1),TrajRef(2),TrajRef(3),'ko')
% plot3(g(1,4),g(2,4),g(4,4),'b*')
% view(70,30)

% -------- Controller --------
if controller == "PD"
psi = 0;

% ------Attitude Controller-----
% Rd = eul2rotm(target_attitude(1:3),'XYZ'); 
psi_d = target_attitude(1:3);
Rz1 = [cos(psi_d(3)) -sin(psi_d(3)) 0;              sin(psi_d(3))  cos(psi_d(3)) 0;            0       0       1      ];
Ry =  [cos(psi_d(2))    0    sin(psi_d(2));             0             1      0;   -sin(psi_d(2))     0   cos(psi_d(2))];
Rz2 = [cos(psi_d(1)) -sin(psi_d(1)) 0;              sin(psi_d(1))  cos(psi_d(1)) 0;            0       0       1      ];
% R0 = eul2rotm(psi','XYZ');
Rd = Rz2 * Ry * Rz1;



R = g(1:3,1:3);

th = pi/4;
alpha= 2; G = eye(3);
rr = r/norm(r);
vv = eul2rotm([pi/2 pi/2 pi/2],'XYZ')*r/norm(r);

e_Ra = 1/2 * unskew(G*Rd'*R - R'*Rd*G);
e_Rb = R'*vv.*rr /(alpha*rr'*R'*vv - alpha*cos(th) );
B =  1 - 1/alpha* log((cos(th) - rr'*R'*vv )/(1 + cos(th)));
A = 1/2* trace(G*(eye(3)-Rd'*R));
err_R =  e_Ra + A*e_Rb;
err_O =  V(1:3)-target_attitude(4:6);
u(1:3) = -0.1*err_R -0.1*err_O;

% ------Trajectory Controller-----
err_pos = TrajRef- r;
err_V = VelRef -V(4:6) ;
u(4:6) =  2*err_V +  0.01*err_pos;





elseif controller == "Backstepping"
psi_d = target_attitude(1:3);
RRef = eul2rotm(psi_d','XYZ');

R = RRef'*R;
r = r -TrajRef;
omega = V(1:3)- target_attitude(4:6);
V = [V(1:3) - target_attitude(4:6); V(4:6)-VelRef];

s = zeros(3,1);
h = s;
for i=1:3
    s = s + a(i)*crossm(R'*e(:,i))*e(:,i);
    h = h + a(i)*crossm(e(:,i))*crossm(omega)*R'*e(:,i); % h = ds/dt
end
l = [s; r];
ldot = [h; R*V(4:6)];
psi = V + K1*l;
u = - ad_conj(V) * II * V -    II * K2 * psi -    II * K1 *  ldot -   II * kappa * [zeros(3,1); R*r];

%Saturation 
% valM=10; valF = 10;
% if u(1) > valM
% u(1)=valM;
% end
% if u(2) > valM
% u(2)=valM;
% end
% if u(3) > valM
% u(3)=valM;
% end
% if u(4) > valF
% u(4)=valF;
% end
% if u(5) > valF
% u(5)=valF;
% end
% if u(6) > valF
% u(6)=valF;
% end

end

end

