%===================================================================================================
%[]SCRIPT NAME: Main_Filter.m
%[]AUTHOR: Matthew Michael Wittal & Gennaro Mangiacapra
%[]CREATED: 10/6/2020
%[]REVISED: 20/6/2020
%===================================================================================================
%[]SCRIPT DESCRIPTION:
%This script simulates pose estimation with UKF for a spacecraft orbiting the irregular moon of 
%Saturn, Pan.
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
Asteroid = "Bennu";
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
N = 1;

dt = 1.2;
%[s]Timestep

ti = 0;
%[s]Initial time

tf = 300;
%[s]Final time

t = ti: dt: tf;
%[s]Timespan

lt = length(t);
%[]Timespan size

th0 = [90, 90, -90]*pi/180; %Pan
% th0 = [-30, -25, 18]*pi/180; %Bennu
R0 = eul2rotm(th0,'XYZ'); 
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
%[km,km/s]Initial satellite state W.R.T. the Asteroid in E.C.E.F. coordina

%% UKF Parameters
measurement_noise_std = [ ones(1,3)* 0.5,...  %Theta [rad]
                          ones(1,3)* 0.1,...  %p [km]
                          ones(1,3)* 0.1*2,...  %w [rad]
                          ones(1,3)* 0.1*2];    %v [km/s]
%[rad,km, rad/s, km/s] Measurement Noise standard deviations

v_k = randn(lt,12)*diag(measurement_noise_std.^2);
%[rad,km, rad/s, km/s] Measurement Noise 

process_noise_std = [ ones(1,3)* 0.005,...  %Theta [rad]
                      ones(1,3)* 0.001 ,...  %p [km]
                      ones(1,3)* 0.001*10,...  %w [rad]
                      ones(1,3)* 0.001*10];    %v [km/s]
%[rad,km, rad/s, km/s] Process Noise standard deviations  

Xbar = eye(4,4);
%Process Noise Mean (always zero)

w_k = randn(lt,12)*diag(process_noise_std.^2);
%[rad,km, rad/s, km/s] Process Noise 

Q = diag(process_noise_std.^2) * eye(12);
cholQ = chol(Q);
%Process noise covariance matrix

R = diag(measurement_noise_std.^2) * eye(12);
%Measurement noise covariance matrix

P0 = eye(12);
%Initial covariance matrix of state estimation error

alpha = [1e-3, 1e-3, 1e-3];
%Sigma point parameters

phi = @SE3_phi;
phi_inv = @SE3_phi_inv;
%Retraction functions

weights = ukf_set_weight(length(P0), length(Q), alpha);
%Weights for the filtering

ukf_state.g = g0; ukf_state.g(1:3,4)=[-2;2;2];
ukf_state.V = [-10;-10;-10; 0.1; 0.1; 0.1];
ukf_P = P0*10;
%UFK initial estimates

ukf_g(:,:,1) =  ukf_state.g; ukf_V(:,:,1) =  ukf_state.V; ukf_P_est(:,:,1) = ukf_P;
%Vectors in which save the estimates

%% ==================================================================================================
%% GUIDANCE SYSTEM
NumPoint = 500;
%Number of points for the discretized orbit

horizon = 5;
%How many points from the closest point on the orbit to reach

e_ref = 0.4;
%[-]Eccentricity

r_p = 50;
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

% figure(101), hold on, grid on
% plot3(target_orbit(:,1),target_orbit(:,2),target_orbit(:,3),'-o')
% plot3(0,0,0,'ok')

hovering =0;
target_orbit(:,1)= -1.5;
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
%% Simulation
%Variables initialization
g_real(:,:,1) = g0;  g_ideal(:,:,1) = g0; 
V_real(:,:,1) = V0; V_ideal(:,:,1) = V0;
Z_real(:, :, 1) = h_sys(g0,V0)+ v_k(1,:)'; 
Z_ideal(:, :, 1) = h_sys(g0,V0);
th(1,1:3)= Z_real(1:3, :, 1);
g_measured(:,:,1) = [expm(crossm(th(1,1:3))),  S_matrix(th(1,1:3))*Z_real(4:6, 1, 1); %expmso3(th)
    zeros(1,3),             1];
V_measured(:,:,1) = Z_real(7:12, :, 1);
    
for i = 2:length(t)

    
    disp("i: "+num2str(i)+"\"+num2str(length(t)))
    % Integration without Noise
    [~, g_, V_, ~] = f_sys(dt, t(i), g_ideal(:,:,i-1), V_ideal(:,:,i-1), zeros(12,1), eye(4,4));
    g_ideal(:,:,i) = g_;
    V_ideal(:,:,i) = V_;
    
    % Integration with Noise
    [~, g_, V_, ~] = f_sys(dt, t(i), g_ideal(:,:,i-1), V_ideal(:,:,i-1) ,w_k(i,1:12)', Xbar);
    g_real(:,:,i) = g_;
    V_real(:,:,i) = V_;
    
    % Measurement Noise
    Z_ideal(:, :, i) = h_sys(g_ideal(:,:,i),V_ideal(:,:,i));
    Z_real(:, :, i) = Z_ideal(:, :, i) + v_k(i,:)';
    
    th(i,1:3)= Z_real(1:3, :, i);
    g_measured(:,:,i) = [expm(crossm(th(i,1:3))),  S_matrix(th(i,1:3))*Z_real(4:6, 1, i); 
           zeros(1,3),             1];
    V_measured(:,:,i) = Z_real(7:12, :, i);
end


for n = 1:lt
    
    
    disp("n: "+num2str(n)+"\"+num2str(lt))
    
    %Propagation Step
    [ukf_state, ukf_P] = ukf_propagation(ukf_state, ukf_P, dt, phi, phi_inv, cholQ, weights, t(n));
    
    %Update Step
    [ukf_state, ukf_P] = ukf_update(ukf_state, ukf_P, Z_real(:,:,n), phi, R, weights);

    %Save Estimates
    ukf_g(:,:,n) = ukf_state.g;
    ukf_V(:,:,n) = ukf_state.V;
    ukf_P_est(:,:,n) = ukf_P;
    std_x(:,n) = sqrt(diag(ukf_P_est(:,:,n)));
    
end

%Store Euler angles from estimation
EulAng_ukf = zeros(lt,3); EulAng_ideal = zeros(lt,3); EulAng_measured = zeros(lt,3);
for i = 1: length(t)
EulAng_ukf(i,:) = rotm2eul(ukf_g(1:3,1:3,i),'XYZ');
EulAng_ideal(i,:) = rotm2eul(g_ideal(1:3,1:3,i),'XYZ');
EulAng_measured(i,:) = rotm2eul(g_measured(1:3,1:3,i),'XYZ');

% asd = vedge_inv(logSE3(ukf_g(:,:,i))); 
% EulAng_ukf(i,:)=asd(1:3);
end

%% Plot States
figure()
subplot(2,3,1), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('x [km]','interpreter','latex');
plot(t, squeeze(g_measured(1, 4, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_g(1, 4, :)),'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(g_ideal(1, 4, :)),'--','color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,2), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('y [km]','interpreter','latex');
plot(t, squeeze(g_measured(2, 4, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_g(2, 4, :)),'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(g_ideal(2, 4, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,3), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('z [km]','interpreter','latex');
plot(t, squeeze(g_measured(3, 4, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_g(3, 4, :)),'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(g_ideal(3, 4, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)
legend('Measured','Filtered','Ideal','interpreter','latex','fontsize',15,'Orientation','horizontal')

subplot(2,3,4), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_z \,[km/s]$','interpreter','latex');
plot(t, squeeze(V_measured(4, :, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, smooth(squeeze(ukf_V(4, :, :))),'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(V_ideal(4, :, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,5), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_y \,[km/s]$','interpreter','latex');
plot(t, squeeze(V_measured(5, :, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(5, :, :)),'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(V_ideal(5, :, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,6), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_z \,[km/s]$','interpreter','latex');
plot(t, squeeze(V_measured(6, :, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(6, :, :)),'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(V_ideal(6, :, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)


set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','ylabel'),'FontSize',36)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])


figure()
% suptitle('Angular Velocity')
subplot(2,3,1), hold on, grid on
ylabel('$\omega_x \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(V_measured(1, :, :))*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(1, :, :))*180/pi,'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(V_ideal(1, :, :))*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,2), hold on, grid on
ylabel('$\omega_y \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(V_measured(2, :, :))*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(2, :, :))*180/pi,'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(V_ideal(2, :, :))*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,3), hold on, grid on
ylabel('$\omega_z \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(V_measured(3, :, :))*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(3, :, :))*180/pi,'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, squeeze(V_ideal(3, :, :)*180/pi),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)
legend('Measured','Filtered','Ideal','interpreter','latex','fontsize',15)

subplot(2,3,4), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_measured(:, 1)*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, EulAng_ukf(:, 1)*180/pi,'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, EulAng_ideal(:, 1)*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,5), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_measured(:, 2)*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, EulAng_ukf(:, 2)*180/pi,'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, EulAng_ideal(:, 2)*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,6), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_measured(:, 3)*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, EulAng_ukf(:, 3)*180/pi,'color',[0.00,0.00,0.00],'linewidth',2)
plot(t, EulAng_ideal(:, 3)*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','ylabel'),'FontSize',36)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])
%% Plot Orbit
S1 = [squeeze(g_ideal(1:3,4,:)); V_ideal(4:6,:)];
S2 = [squeeze(ukf_g(1:3,4,:)); ukf_V(4:6,:)];
S3 = [squeeze(g_measured(1:3,4,:)); V_measured(4:6,:)];

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


%% UKF Perf
figure()
subplot(2,3,1), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$\epsilon_x$ [km]','interpreter','latex');
plot(t, squeeze(ukf_g(1, 4, :))-squeeze(g_ideal(1, 4, :)),'linewidth',2)
plot(t,3*std_x(4,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(4,:),'--','linewidth',2);

subplot(2,3,2), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$\epsilon_y$ [km]','interpreter','latex');
plot(t, squeeze(ukf_g(2, 4, :))-squeeze(g_ideal(2, 4, :)),'linewidth',2)
plot(t,3*std_x(5,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(5,:),'--','linewidth',2);

subplot(2,3,3), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$\epsilon_z$ [km]','interpreter','latex');
h1=plot(t, squeeze(ukf_g(3, 4, :))-squeeze(g_ideal(3, 4, :)),'linewidth',2)
h2=plot(t,3*std_x(6,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
h3=plot(t,-3*std_x(6,:),'--','linewidth',2);
legend([h2 h3],'$+3\sigma$ bound','$-3\sigma$ bound','interpreter','latex','Fontsize',12,'location','southeast','Orientation','horizontal')

subplot(2,3,4), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$\epsilon_{v_x}$ [km/s]','interpreter','latex');
plot(t, smooth(squeeze(ukf_V(4, :, :)))-squeeze(V_ideal(4, :, :)),'linewidth',2)
plot(t,3*smooth(std_x(10,:)),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*smooth(std_x(10,:)),'--','linewidth',2);

subplot(2,3,5), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$\epsilon_{v_y}$ [km/s]','interpreter','latex');
plot(t, squeeze(ukf_V(5, :, :))-squeeze(V_ideal(5, :, :)),'linewidth',2)
plot(t,3*std_x(11,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(11,:),'--','linewidth',2);

subplot(2,3,6), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$\epsilon_{v_z}$ [km/s]','interpreter','latex');
plot(t, squeeze(ukf_V(6, :, :))-squeeze(V_ideal(6, :, :)),'linewidth',2)
plot(t,3*std_x(12,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(12,:),'--','linewidth',2);


set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','ylabel'),'FontSize',36)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])


figure()
subplot(2,3,1), hold on, grid on
ylabel('$\epsilon_{\omega_x}$ [deg/s]','interpreter','latex');
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(ukf_V(1, :, :))*180/pi-squeeze(V_ideal(1, :, :))*180/pi,'linewidth',2)
plot(t,3*std_x(7,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(7,:)*180/pi,'--','linewidth',2);

subplot(2,3,2), hold on, grid on
ylabel('$\epsilon_{\omega_y}$ [deg/s]','interpreter','latex');
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(ukf_V(2, :, :))*180/pi-squeeze(V_ideal(2, :, :))*180/pi,'linewidth',2)
plot(t,3*std_x(8,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(8,:)*180/pi,'--','linewidth',2);


subplot(2,3,3), hold on, grid on
ylabel('$\epsilon_{\omega_z}$ [deg/s]','interpreter','latex');
xlabel('t [s]','interpreter','latex');
h1=plot(t, squeeze(ukf_V(3, :, :))*180/pi-squeeze(V_ideal(3, :, :))*180/pi,'linewidth',2)
h2=plot(t,3*std_x(9,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
h3=plot(t,-3*std_x(9,:)*180/pi,'--','linewidth',2);
legend([h2 h3],'$+3\sigma$ bound','$-3\sigma$ bound','interpreter','latex','Fontsize',12,'location','southeast')

subplot(2,3,4), hold on, grid on
ylabel('$\epsilon_{\Theta_1}$ [deg/s]','interpreter','latex');
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_ukf(:, 1)*180/pi- EulAng_ideal(:, 1)*180/pi,'linewidth',2)
plot(t,3*std_x(1,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(1,:)*180/pi,'--','linewidth',2);

subplot(2,3,5), hold on, grid on
ylabel('$\epsilon_{\Theta_2}$ [deg/s]','interpreter','latex');
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_ukf(:, 2)*180/pi- EulAng_ideal(:, 2)*180/pi,'linewidth',2)
plot(t,3*std_x(2,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(2,:)*180/pi,'--','linewidth',2);

subplot(2,3,6), hold on, grid on
ylabel('$\epsilon_{\Theta_3}$ [deg/s]','interpreter','latex');
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_ukf(:, 3)*180/pi- EulAng_ideal(:, 3)*180/pi,'linewidth',2)
plot(t,3*std_x(3,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(3,:)*180/pi,'--','linewidth',2);

set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','ylabel'),'FontSize',36)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])

%% Control Effort
figure(), 

subplot(2,1,1),grid on, hold on,
plot(t,UU(1,1:204:end-1),'linewidth',1)
plot(t,UU(2,1:204:end-1),'linewidth',1)
plot(t,UU(3,1:204:end-1),'linewidth',1)
ylabel('$u_M \, [kN m]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
legend('$u_{M_x}$','$u_{M_y}$','$u_{M_z}$','interpreter','latex','fontsize',15)

subplot(2,1,2),grid on, hold on,
plot(t,UU(4,1:204:end-1),'linewidth',1)
plot(t,UU(5,1:204:end-1),'linewidth',1)
plot(t,UU(6,1:204:end-1),'linewidth',1)
ylabel('$u_F \, [kN]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
legend('$u_{F_x}$','$u_{F_y}$','$u_{F_z}$','interpreter','latex','fontsize',15)

set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])

return
%% Video Visualization
BB = [ukf_g(1, 4, :) ukf_g(2, 4, :) ukf_g(3, 4, :)]; % Position data
A = ukf_g(1:3, 1:3, :) ;                     % Orientation data

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


% ukf_V(4, :, :) = smooth(ukf_V(4, :, :),10);
% ukf_V(5, :, :) = smooth(ukf_V(5, :, :),10);
% ukf_V(6, :, :) = smooth(ukf_V(6, :, :),10);
% std_x(10,:) =  smooth(std_x(10,:) ,100);
% std_x(11,:)=  smooth(std_x(11,:) ,100);
% std_x(12,:)=  smooth(std_x(12,:) ,100);
