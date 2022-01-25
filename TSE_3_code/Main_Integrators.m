%===================================================================================================
%[]SCRIPT NAME: Main_Integrators.m
%[]AUTHOR: Matthew Michael Wittal & Gennaro Mangiacapra
%[]CREATED: 10/6/2020
%[]REVISED: 23/7/2020
%===================================================================================================
%[]SCRIPT DESCRIPTION:
%This script simulates TSE3 and point mass numerical integration schemes 
%for a spacecraft orbiting the irregular asteroid.
%===================================================================================================
%[]COMMENTS:
%'WRT' stands for "with respect to". 'PCI' stands for "Asteroid-Centered Inertial". 'PCPF' stands for
%"Asteroid-Centered Asteroid-Fixed". 'DCM' stands for "Direction Cosine Matrix".
%===================================================================================================
tic;
%[]Starts program timer.

clc;
%[]Clears command window.

clear;
%[]Clears variable workspace.

format long g;
%[]Adjusts the format of the command window output.

format compact;
%[]Adjusts the format of the command window output.

close all;
%[]Closes all windows.

%%==================================================================================================
%% SPACECRAFT AND ORBITAL PROPERTIES
global  II lt J m control mu Asteroid C20 C22 Wp

J = diag([4.97 6.16 8.37])/1e6;
%[kg.m^2]Inertia tensor.

m = 60;
%[kg]Mass of the spacecraft.

II = [J         zeros(3);
    zeros(3)  m*eye(3)];
%[kg.m^2,~,~,kg]Moment of Inertia tensor.

control = 'no';
%Enable/Disable control

%%==================================================================================================
%% ASTEROID SELECTION
Asteroid = "Pan";
%[-]Define Asteroid (Pan, Beenu, Itokawa)

if Asteroid == "Pan"
    
    mu = 4.95*1e15 *  6.67408*1e-11 *1e-9;
    %[kg/km^3]gravitational parameter of Asteroid.
    
    Wp = 2 * pi / 49684.3812 * [0; 0; 1];
    %[rad/s]Angular velocity of Pan in P.C.E.F. coordinates.
    
    l1 = 34.4;
    l2 = 31.4;
    l3 = 20.8;
    %[km]Semimajor axes of oblate spheroid Pan.
    
    r0 = 75*[1.1 ; 0; 0.08];
    %[km]Initial position.

elseif Asteroid == "Bennu"
    
    mu = 5.2*1e-9;
    %[km^3/s^2]gravitational parameter of Bemmu.
    
    Wp = 2 * pi / (4.297*60*60) * [0; 0; 1];
    %[rad/s]Angular velocity of Bennu in P.C.E.F. coordinates.
    
    l1 = 0.565;
    l2 = 0.535;
    l3 = 0.508;
    %[km]Semimajor axes of oblate spheroid Itokawa.
    
    r0 = 2.2*[sqrt(11/12) ; 0; sqrt(1/12)];
    %[km]Initial position.
    
elseif Asteroid == "Itokawa"
    
    mu = 3.51e10*6.67e-11*1e-9;
    %[km^3/s^2]gravitational parameter of Itokawa.
    
    Wp  = 2 * pi / (12.132*60*60) * [0; 0; 1];
    %[rad/s]Angular velocity of Itokawa in P.C.E.F. coordinates.
    
    l1 = 0.535;
    l2 = 0.294;
    l3 = 0.210;
    %[km]Semimajor axes of oblate spheroid Itokawa.

    r0 = 1.5*[sqrt(11/12) ; 0; sqrt(1/12)];
    %[km]Initial position.
    
end

ro = l1/2;
alpha = l1/l1;
beta = l2/l1;
gamma = l3/l1;
%[]Dimensionless axes.

C20 = 1/(5*ro^2)*(gamma^2 -(alpha^2+beta^2)/2);
C22 = 1/(20*ro^2)*(alpha^2-beta^2);
%[]Harmonic Terms
 
%% TIMESPAN AND INITIAL CONDITIONS
dt = 3600;
%[s]Timestep

ti = 0;
%[s]Initial time

tf = 86400*20;
%[s]Final time

t = ti: dt: tf;
%[s]Timespan

lt = length(t);
%[]Timespan size

th0 = [0, 0, 0]*pi/180;
R0 = eul2rotm(th0,'XYZ'); %
%[rad]Initial attitude DCM.

w0 = [0.000001;  0.000001;  0.000001]*pi/180; 
%[rad/s]Initial attitude rate.

% v0 =  1 * sqrt(mu/norm(r0))*cross(r0,[0;0;1])/norm(cross(r0,[0,0,1]));
v0 =  1 * sqrt(mu/norm(r0))*[0.4;1;0];
%[km/s]Initial velocity.

g0 = [R0 r0; zeros(1,3) 1];
%[rad,km,~,~]Initial Pose in SE(3).

V0 = [w0; v0];
%[rad/s,km/s]Initial pose rate.

X0 = vedge_inv( logSE3( g0 ) );
%[rad,km] X0 corresponds to g0 in the Lie Algebra se(3):

So = [r0; v0];
%[km,km/s]Initial satellite state W.R.T. the Asteroid in E.C.E.F. coordinates.

%%==================================================================================================
%% INTEGRATIONS
% Integrator on TSE(3)
[psi, g, V, tet] = rk4_dynamicsTSE3(dt, ti, tf, g0, V0);

% Point mass integrator
options = odeset('RelTol',10^-12, 'AbsTol', 10^-12);
[t3,S3] = ode45(@(t2,S2)PointMassModel_Pan(t2,S2),t,So,options);


%Save Cartesian state
S1 = [squeeze(g(1:3,4,:)); V(4:6,:)];
S3=S3';

%Compute gravitational forces[kN] and moments [MNm] in the Asteroid reference frame
F_SE3 = zeros(lt,3); M_SE3 = zeros(lt,3); 
F_se3 = zeros(lt,3); M_se3 = zeros(lt,3); 
F_Euclid = zeros(lt,3); M_Euclid = zeros(lt,3);

for i = 1: length(t)
[F,M] = Pan_GravityModel(g(1:3,4,i),g(1:3,1:3,i),V(:,i),t(i));   
F_SE3(i,:) = F';
M_SE3(i,:) = M';
 

[F,M] = Pan_GravityModel(S3(1:3,i),eye(3,3),[zeros(3,1);S3(4:6,i)],t(i));  
F_Euclid(i,:) = F';
M_Euclid(i,:) = M';

end

%Compute Euler Angles for SE(3) model
EulAng = zeros(lt,3);
for i = 1: length(t)
EulAng(i,:) = rotm2eul(g(1:3,1:3,i),'XYZ');
end
%%==================================================================================================
%% Plot 3D Orbit and Forces/Moments
Extent = 5;
%[Asteroid radii]Axes extent.

Projection= 0;
%[-]Enable/Disable projection plot 1/0 (ON/OFF), works only with two orbits

AXES= 0;
%[-]Enable/Disable asteroid axes 1/0 (ON/OFF)

PlotOrbit(AXES,Extent,Asteroid,Projection,S1(1:3,:),'Rigid Body Motion', ...
       S3(1:3,:),'Point Body Motion')
%[]Plots the orbit in three (3) dimensions.

figure()
subplot(2,3,1),hold on, grid on
plot(t,F_SE3(:,1),'k-','linewidth',1)
plot(t,F_Euclid(:,1),'b-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$F_x \, [kN]$','interpreter','latex');

subplot(2,3,2),hold on, grid on
plot(t,F_SE3(:,2),'k-','linewidth',1)
plot(t,F_Euclid(:,2),'b-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$F_y \, [kN]$','interpreter','latex');

subplot(2,3,3),hold on, grid on
plot(t,F_SE3(:,3),'k-','linewidth',1)
plot(t,F_Euclid(:,3),'b-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$F_y \, [kN]$','interpreter','latex');
legend('Rigid Body','Point Mass','interpreter','latex','fontsize',15)

subplot(2,3,4),hold on, grid on
plot(t,M_SE3(:,1),'k-','linewidth',1)
plot(t,M_Euclid(:,1),'b-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$M_x \, [kNm]$','interpreter','latex');

subplot(2,3,5),hold on, grid on
plot(t,M_SE3(:,2),'k-','linewidth',1)
plot(t,M_Euclid(:,2),'b-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$M_y \, [kNm]$','interpreter','latex');

subplot(2,3,6),hold on, grid on
plot(t,M_SE3(:,3),'k-','linewidth',1)
plot(t,M_Euclid(:,3),'b-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$M_z \, [kNm]$','interpreter','latex');

set(findall(gcf,'-property','FontSize'),'FontSize',17)
%[kN] Plots gravitational force components for the three integration schemes
%[MNm] Plots gravitational moments components for the three integration schemes

%% Plot States
figure(), 
subplot(2,2,1),hold on, grid on
plot(t, squeeze(g(1:3, 4, :)),'linewidth',2)
xlabel('t [s]','interpreter','latex');
ylabel('Position [km]','interpreter','latex');
legend('x','y','z','interpreter','latex')
set(gca,'FontSize',17)

subplot(2,2,2),hold on, grid on
plot(t, squeeze(V(1:3, :, :))*180/pi,'linewidth',2)
xlabel('t [s]','interpreter','latex');
ylabel('Angular Velocity [deg/s]','interpreter','latex');
legend('$\omega_x$','$\omega_y$','$\omega_z$','interpreter','latex')
set(gca,'FontSize',17)

subplot(2,2,3),hold on, grid on
plot(t, squeeze(V(4:6, :, :)),'linewidth',2)
xlabel('t [s]','interpreter','latex');
ylabel('Linear Velocity [km/s]','interpreter','latex');
legend('$v_x$','$v_y$','$v_z$','interpreter','latex')
set(gca,'FontSize',17)

subplot(2,2,4),hold on, grid on
plot(t,EulAng(:,1)*180/pi,'linewidth',2)
plot(t,EulAng(:,2)*180/pi,'linewidth',2), 
plot(t,EulAng(:,3)*180/pi,'linewidth',2), 
xlabel('t [s]','interpreter','latex');
ylabel('Euler angles [deg]','interpreter','latex');
legend('$\Theta_1$','$\Theta_2$','$\Theta_3$','interpreter','latex','fontsize',15)
set(gca,'FontSize',17)


%% Kepler Elements Plot
for i = 1:lt
[a_d,e_d,i_d,Omega_d,omega_d,M_d,E_d,theta_d,l_d,Pi_d,u_d]= cart2kepl(S1(1,i),S1(2,i),S1(3,i),S1(4,i),S1(5,i),S1(6,i),mu);

a_m1(i)=a_d;
ecc1(i)=e_d;
inc1(i)=i_d*pi/180;
w1(i)= omega_d*pi/180;
AN1(i)= Omega_d*pi/180;
th1(i)= l_d*pi/180;
end

for i = 1:lt
[a_d,e_d,i_d,Omega_d,omega_d,M_d,E_d,theta_d,l_d,Pi_d,u_d]= cart2kepl(S3(1,i),S3(2,i),S3(3,i),S3(4,i),S3(5,i),S3(6,i),mu);

a_m3(i)=a_d;
ecc3(i)=e_d;
inc3(i)=i_d*pi/180;
w3(i)= omega_d*pi/180;
AN3(i)= Omega_d*pi/180;
th3(i)= l_d*pi/180;
end

figure()
subplot(3,2,1), hold on, grid on
plot(t,a_m1,'linewidth',1)
plot(t,a_m3,'linewidth',1)
xlabel('$t \ [s]$','interpreter','latex','fontsize',13)
ylabel('$a \ [m]$','interpreter','latex','fontsize',13)
title('Semi major axis','interpreter','latex','fontsize',13)

subplot(3,2,2), hold on, grid on
plot(t,ecc1,'linewidth',1)
plot(t,ecc3,'linewidth',1)
xlabel('$t \ [s]$','interpreter','latex','fontsize',13)
ylabel('$e \ [-]$','interpreter','latex','fontsize',13)
title('Eccentricity','interpreter','latex','fontsize',13)
legend('Rigid Body','Point Mass')

subplot(3,2,3), hold on, grid on
plot(t,(inc1*180/pi),'linewidth',1)
plot(t,(inc3*180/pi),'linewidth',1)
xlabel('$t \ [s]$','interpreter','latex','fontsize',13)
ylabel('$i \ [deg]$','interpreter','latex','fontsize',13)
title('Inclination','interpreter','latex','fontsize',13)

subplot(3,2,4), hold on, grid on
plot(t,(w1),'linewidth',1)
plot(t,(w3),'linewidth',1)
xlabel('$t \ [s]$','interpreter','latex','fontsize',13)
ylabel('$\omega \ [rad]$','interpreter','latex','fontsize',13)
title('Argument of periapsis','interpreter','latex','fontsize',13)

subplot(3,2,5), hold on, grid on
plot(t,(AN1),'linewidth',1)
plot(t,(AN3),'linewidth',1)
xlabel('$t \ [s]$','interpreter','latex','fontsize',13)
ylabel('$\Omega \ [rad]$','interpreter','latex','fontsize',13)
title('Longitude of ascending node','interpreter','latex','fontsize',13)

subplot(3,2,6), hold on, grid on
plot(t,(th1*180/pi),'linewidth',1)
plot(t,(th3*180/pi),'linewidth',1)
xlabel('$t \ [s]$','interpreter','latex','fontsize',13)
ylabel('$\theta \ [rad]$','interpreter','latex','fontsize',13)
title('True anomaly','interpreter','latex','fontsize',13)

return
%% Video Visualization
BB = [g(1, 4, :) g(2, 4, :) g(3, 4, :)]; % Position data
A = g(1:3, 1:3, :) ;                     % Orientation data

myvizualization_mod(t,BB,A)

%% Ground Track
PlotGroundTrack_New(t,S1(1:3,:),S3(1:3,:),Asteroid)
AnimateGroundTrack_New(t,S1(1:3,:),S3(1:3,:),dt,Asteroid)
%[]Ground track plot and animation