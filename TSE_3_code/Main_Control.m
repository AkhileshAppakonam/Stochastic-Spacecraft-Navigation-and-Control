%===================================================================================================
%[]AUTHOR: Matthew Michael Wittal & Gennaro Mangiacapra
%[]CREATED: 10/6/2020
%[]REVISED: 23/7/2020
%===================================================================================================
%[]SCRIPT DESCRIPTION:
%This script simulates the spacecraft motion around the asteroid. A
%backstepping control logic is employed to control both the attitude and
%the trajectory
%===================================================================================================
tic;
%[]Starts program timer.

clc;
%[]Clears command window.

clear all;
clear global;
%[]Clears variable workspace.

format long g;
%[]Adjusts the format of the command window output.

format compact;
%[]Adjusts the format of the command window output.

close all;
%[]Closes all windows.

%%==================================================================================================
%% SPACECRAFT AND ORBITAL PROPERTIES
global  II lt K1 K2 kappa e a J m UU control mu target_orbit controller target_attitude Asteroid horizon C20 C22

J = diag([4.97 6.16 8.37])/1e6;  % kg.m^2
%[kg.m^2]Inertia tensor.

m = 60;
%[kg]Mass of the spacecraft.

II = [J         zeros(3);
    zeros(3)  m*eye(3)];
%[kg.m^2,~,~,kg]Moment of Inertia tensor.


UU = zeros(6,1);
%Control effort variable

control = 'yes';
%Enable/Disable control (yes/no)

controller = 'Backstepping';
%Type of controller (PD/Backstepping)

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

elseif Asteroid == "Bennu"
    
    mu = 5.2*1e-9;
    %[km^3/s^2]gravitational parameter of Bemmu.
    
    Wp = 2 * pi / (4.297*60*60) * [0; 0; 1];
    %[rad/s]Angular velocity of Bennu in P.C.E.F. coordinates.
    
    l1 = 0.565;
    l2 = 0.535;
    l3 = 0.508;
    %[km]Semimajor axes of oblate spheroid Itokawa.
    
elseif Asteroid == "Itokawa"
    
    mu = 3.51e10*6.67e-11*1e-9;
    %[km^3/s^2]gravitational parameter of Itokawa.
    
    Wp  = 2 * pi / (12.132*60*60) * [0; 0; 1];
    %[rad/s]Angular velocity of Itokawa in P.C.E.F. coordinates.
    
    l1 = 0.535;
    l2 = 0.294;
    l3 = 0.210;
    %[km]Semimajor axes of oblate spheroid Itokawa.

end

ro = l1/2;
alpha = l1/l1;
beta = l2/l1;
gamma = l3/l1;
%[]Dimensionless axes.

C20 = 1/(5*ro^2)*(gamma^2 -(alpha^2+beta^2)/2);
C22 = 1/(20*ro^2)*(alpha^2-beta^2);
%[]Harmonic Terms
%%==================================================================================================
%% TIMESPAN AND INITIAL CONDITIONS
dt = 1.2;
%[s]Timestep

ti = 0;
%[s]Initial time

tf = 2000;
%[s]Final time

t = ti: dt: tf;
%[s]Timespan

lt = length(t);
%[]Timespan size

th0 = [90, 90, -90]*pi/180; %Pan
% th0 = [-30, -25, 18]*pi/180; %Bennu
R0 = eul2rotm(th0,'XYZ'); %
%[rad]Initial attitude DCM.

r0 = [65; 10; 15]; %Pan
% r0 = [1.6; -1.3; -1.0]; %Bennu
%[km]Initial position.

w0 = [15;  3;  9]*pi/180; %Pan
% w0 = [120;  30;  90]*pi/180;  %Bennu
%[rad/s]Initial attitude rate.

v0 =  1 * sqrt(mu/norm(r0))*[0.4;1;0];  %Pan
% v0 =  1 * sqrt(mu/norm(r0))*[0.78;1;0.3];  %Bennu
%[km/s]Initial velocity.

%Backstepping Controller Parameters:
% kappa = 0.0002; %Pan
kappa = 0.002; %Bennu
k11 = 0.2; 
k12 = 0.1; %Pan
% k12 = 1.0;  %Bennu
k21 = 1.0; 
k22 = 0.01;%Pan
k22 = 0.01; %Bennu
K1 = blkdiag(k11*eye(3), k12*eye(3));
K2 = blkdiag(k21*eye(3), k22*eye(3));
a = [1.2; 1.1; 1]; % a1 > a2 > a3 >= 1
e = eye(3);

g0 = [R0 r0; zeros(1,3) 1];
%[rad,km,~,~]Initial Pose in SE(3).

V0 = [w0; v0];
%[rad/s,km/s]Initial pose rate.

X0 = vedge_inv( logSE3( g0 ) );
%[rad,km] X0 corresponds to g0 in the Lie Algebra se(3):

So = [r0; v0];
%[km,km/s]Initial satellite state W.R.T. the Asteroid in E.C.E.F. coordinates.

%% ==================================================================================================
%% GUIDANCE SYSTEM
NumPoint = 500;
%Number of points for the discretized orbit

horizon = 5;
%How many points from the closest point on the orbit to reach

e_ref = 1e-3;
%[-]Eccentricity

r_p = 60;
%[km]Distance at periapsis

a_ref = r_p*(1+e_ref)/(1-e_ref^2);
%[km]Semi-major axis

i_ref = 25;
%[deg]Inclination

Omega_ref = 0;
%[deg]RAAN

omega_ref = 0;
%[deg]Argument of perigee

n_ref = sqrt(mu/a_ref^3);
%[rad/s]Mean motion

l_ref = 0;
%[deg]True longitude

Pi_ref = 0;
%[deg]Longitude of perigee 

u_ref = 0;
%[deg]Argument of latitude

tau=0;
%[s]Time passage to perigee

tt = linspace(0,round(2*pi/n_ref)*180/pi,NumPoint);
%[s]Time window on which create the orbit dataset

%Preallocation
i=1;
target_orbit = zeros(length(tt),6);
target_attitude = zeros(6,1);

for k = tt

M_ref = n_ref* (k-tau);
%[deg]Mean anomaly

[x,y,z,xdot,ydot,zdot] = kepl2cart(a_ref,e_ref,i_ref,Omega_ref,omega_ref,M_ref,l_ref,Pi_ref,u_ref,mu);

%Store Cartesian state
target_orbit(i,1:3)  = [x,y,z];
target_orbit(i,4:6)  = [xdot,ydot,zdot];

i=i+1;
end

hovering=0;
% target_orbit(:,1)= -1.5;
% target_orbit(:,2)= 1;
% target_orbit(:,3)= 0.5;
% target_orbit(:,4)= 0;
% target_orbit(:,5)= 0;
% target_orbit(:,6)= 0;
% hovering=1;
% horizon = horizon*0;
%Uncomment for Hovering problem

target_attitude(1:3,1)= [12;5;9]*pi/180; %Pan
% target_attitude(1:3,1)= [30;20;10]*pi/180; %Bennu
%[rad]Target Attitude

target_attitude(4:6,1)= [0.0;0.0;0.0];
%[rad/s]Target Angular Velocity 
%%==================================================================================================
%% INTEGRATIONS
%Preallocation
g = zeros(4,4,lt); V = zeros(6,1,lt); Z = zeros(12,1,lt);
g(:,:,1) = g0; V(:,:,1) = V0; Z(:,:,1) = [X0;V0];
for i = 2:lt

    %Attitude Doublet
%     if t(i)< 800
%         target_attitude(1:3,1)= [12;12;12]*pi/180; %Attitude
%     elseif 800<=t(i) && t(i)<2000
%         target_attitude(1:3,1)= -[7;7;7]*pi/180; %Attitude
%     elseif t(i)>=2000
%         target_attitude(1:3,1)= [9;9;9]*pi/180; %Attitude
%     end


    %Show current index
    disp("i: "+num2str(i)+"\"+num2str(length(t)))

    %Integration
    [~, g_, V_, ~] = f_sys(dt, t(i), g(:,:,i-1), V(:,:,i-1), zeros(12,1), eye(4,4));
    g(:,:,i) = g_;
    V(:,:,i) = V_;

    %Measurement
    Z(:, :, i) = h_sys(g(:,:,i),V(:,:,i));
    
end

%Save Cartesian state
S1 = [squeeze(g(1:3,4,:)); V(4:6,:)];

%Compute gravitational forces[kN] and moments [kNm] in the Pan reference frame
F_SE3 = zeros(lt,3); M_SE3 = zeros(lt,3); 
for i = 1: length(t)
[F,M] = Pan_GravityModel(g(1:3,4,i),g(1:3,1:3,i),V(:,i),t(i));   
F_SE3(i,:) = F';
M_SE3(i,:) = M';
end

%Store Euler angles from pose
EulAng = zeros(lt,3);
for i = 1: length(t)
EulAng(i,:) = rotm2eul(g(1:3,1:3,i),'XYZ');
end

%% Plot 3D Orbit and Forces/Moments
Extent = 5;
%[Asteroid radii]Axes extent.

AXES= 0;
%[-]Enable/Disable asteroid axes 1/0 (ON/OFF)

Projection= 1;
%[-]Enable/Disable projection plot 1/0 (ON/OFF)
if hovering==0
PlotOrbit(AXES,Extent,Asteroid,Projection,S1(1:3,:),'Spacecraft Orbit',target_orbit(:,1:3)','Target Orbit');
else
PlotOrbit(AXES,Extent,Asteroid,Projection,S1(1:3,:),'Spacecraft Orbit',target_orbit(:,1:3)','Target Position');
end
%[-]Plots the orbit in three (3) dimensions.

figure()
subplot(2,3,1),hold on, grid on
plot(t,F_SE3(:,1),'k-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$F_x \, [kN]$','interpreter','latex');

subplot(2,3,2),hold on, grid on
plot(t,F_SE3(:,2),'k-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$F_y \, [kN]$','interpreter','latex');

subplot(2,3,3),hold on, grid on
plot(t,F_SE3(:,3),'k-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$F_y \, [kN]$','interpreter','latex');

subplot(2,3,4),hold on, grid on
plot(t,M_SE3(:,1),'k-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$M_x \, [kNm]$','interpreter','latex');

subplot(2,3,5),hold on, grid on
plot(t,M_SE3(:,2),'k-','linewidth',1)
xlabel('t [s]','interpreter','latex');
ylabel('$M_y \, [kNm]$','interpreter','latex');

subplot(2,3,6),hold on, grid on
plot(t,M_SE3(:,3),'k-','linewidth',1)
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

subplot(2,2,2),hold on, grid on
plot(t, squeeze(V(1:3, :, :))*180/pi,'linewidth',2)
xlabel('t [s]','interpreter','latex');
ylabel('Angular Velocity [deg/s]','interpreter','latex');
legend('$\omega_x$','$\omega_y$','$\omega_z$','interpreter','latex')

subplot(2,2,3),hold on, grid on
plot(t, squeeze(V(4:6, :, :)),'linewidth',2)
xlabel('t [s]','interpreter','latex');
ylabel('Linear Velocity [km/s]','interpreter','latex');
legend('$v_x$','$v_y$','$v_z$','interpreter','latex')

subplot(2,2,4),hold on, grid on
plot(t,EulAng(:,1)*180/pi,'linewidth',2)
plot(t,EulAng(:,2)*180/pi,'linewidth',2), 
plot(t,EulAng(:,3)*180/pi,'linewidth',2), 
xlabel('t [s]','interpreter','latex');
ylabel('Euler angles [deg]','interpreter','latex');
legend('$\Theta_1$','$\Theta_2$','$\Theta_3$','interpreter','latex','fontsize',15)

set(findall(gcf,'-property','FontSize'),'FontSize',17)
%% Control Effort
figure()
subplot(2,1,1),grid on, hold on,
plot(t,UU(1,1:4:end),'linewidth',1)
plot(t,UU(2,1:4:end),'linewidth',1)
plot(t,UU(3,1:4:end),'linewidth',1)
ylabel('$u_M \, [kN m]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
legend('$u_{M_x}$','$u_{M_y}$','$u_{M_z}$','interpreter','latex','fontsize',15)

subplot(2,1,2),grid on, hold on,
plot(t,UU(4,1:4:end),'linewidth',1)
plot(t,UU(5,1:4:end),'linewidth',1)
plot(t,UU(6,1:4:end),'linewidth',1)
ylabel('$u_F \, [kN]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
legend('$u_{F_x}$','$u_{F_y}$','$u_{F_z}$','interpreter','latex','fontsize',15)

set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])

return
%% Video Visualization
BB = [g(1, 4, :) g(2, 4, :) g(3, 4, :)]; % Position data
A = g(1:3, 1:3, :) ;                     % Orientation data

myvizualization_mod(t,BB,A)

%% Ground Track
PlotGroundTrack_New(t,S1(1:3,:),S2(1:3,:),Asteroid)
AnimateGroundTrack_New(t,S1(1:3,:),S2(1:3,:),dt,Asteroid)
%[]Ground track plot and animation

%     %Attitude Doublet
%     if t(i)< 800
%         target_attitude(1:3,1)= [12;12;12]*pi/180; %Attitude
%     elseif 800<=t(i) && t(i)<2000
%         target_attitude(1:3,1)= -[7;7;7]*pi/180; %Attitude
%     elseif t(i)>=2000
%         target_attitude(1:3,1)= [9;9;9]*pi/180; %Attitude
%     end

