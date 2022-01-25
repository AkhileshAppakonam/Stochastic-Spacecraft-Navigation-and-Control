%===================================================================================================
%[]AUTHOR: Matthew Michael Wittal & Gennaro Mangiacapra
%[]CREATED: 10/6/2020
%[]REVISED: 23/7/2020
%===================================================================================================
%[]SCRIPT DESCRIPTION:
%This script simulates pose estimation with UKF for a spacecraft orbiting 
% the irregular moon asteroid
%===================================================================================================

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
global  II lt N J m control mu Asteroid C20 C22 Wp

J = diag([4.97 6.16 8.37]);
%[kg.m^2]Inertia tensor.

m = 60;
%[kg]Mass of the spacecraft.

II = [J         zeros(3);
    zeros(3)  m*eye(3)];
%[kg.m^2,~,~,kg]Moment of Inertia tensor.

control = 'no';
%Enable/Disable control

mu = 4.95*1e15 *  6.67408*1e-11 *1e-9;
%[kg/km^3]gravitational parameter of Asteroid.

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

%%==================================================================================================
%% TIMESPAN AND INITIAL CONDITIONS
N = 1;

dt = 500;
%[s]Timestep

ti = 0;
%[s]Initial time

tf = 86400*1;
%[s]Final time

t = ti: dt: tf;
%[s]Timespan

lt = length(t);
%[]Timespan size

th0 = [0, 0, 0]*pi/180;
R0 = eul2rotm(th0,'XYZ'); %
%[rad]Initial attitude DCM.

r0 = [65; 0; 6.5];
%[km]Initial position.

w0 = [0.000001;  0.000001;  0.000001]*pi/180; 
%[rad/s]Initial attitude rate.

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

%% UKF Parameters
measurement_noise_std = [ ones(1,3)* 0.1,...  %Theta [rad]
                          ones(1,3)* 0.1,...  %p [km]
                          ones(1,3)* 0.1,...  %w [rad]
                          ones(1,3)* 0.1];    %v [km/s]
%[rad,km, rad/s, km/s] Measurement Noise standard deviations

v_k = randn(lt,12)*diag(measurement_noise_std.^2);
%[rad,km, rad/s, km/s] Measurement Noise 

process_noise_std = [ ones(1,3)* 0.001,...  %Theta [rad]
                      ones(1,3)* 0.001 ,...  %p [km]
                      ones(1,3)* 0.001,...  %w [rad]
                      ones(1,3)* 0.001];    %v [km/s]
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

ukf_state.g = g0; ukf_state.V = V0; ukf_P = P0;
%UFK initial estimates

ukf_g(:,:,1) =  ukf_state.g; ukf_V(:,:,1) =  ukf_state.V; ukf_P_est(:,:,1) = ukf_P;
%Vectors in which save the estimates

%%==================================================================================================
%% Simulation 
%Variables initialization
g_real(:,:,1) = g0;  g_ideal(:,:,1) = g0; 
V_real(:,:,1) = V0; V_ideal(:,:,1) = V0;
Z_real(:, :, 1) = h_sys(g0,V0)+ v_k(1,:)'; 
Z_ideal(:, :, 1) = h_sys(g0,V0);
th(1,1:3)= Z_ideal(1:3, :, 1);
g_measured(:,:,1) = [expm(crossm(th(1,1:3))),  S_matrix(th(1,1:3))*Z_ideal(4:6, 1, 1); %expmso3(th)
    zeros(1,3),             1];
V_measured(:,:,1) = Z_ideal(7:12, :, 1);

for i = 2:length(t)

    %Show current index
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

    %Show current index
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
plot(t, squeeze(ukf_g(1, 4, :)),'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(g_ideal(1, 4, :)),'--','color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,2), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('y [km]','interpreter','latex');
plot(t, squeeze(g_measured(2, 4, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_g(2, 4, :)),'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(g_ideal(2, 4, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,3), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('z [km]','interpreter','latex');
plot(t, squeeze(g_measured(3, 4, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_g(3, 4, :)),'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(g_ideal(3, 4, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)
legend('Measured','Filtered','Ideal','interpreter','latex','fontsize',15)

subplot(2,3,4), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_z \,[km/s]$','interpreter','latex');
plot(t, squeeze(V_measured(4, :, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(4, :, :)),'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(V_ideal(4, :, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,5), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_y \,[km/s]$','interpreter','latex');
plot(t, squeeze(V_measured(5, :, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(5, :, :)),'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(V_ideal(5, :, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,6), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_z \,[km/s]$','interpreter','latex');
plot(t, squeeze(V_measured(6, :, :)),'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(6, :, :)),'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(V_ideal(6, :, :)),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)


set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])



figure()
% suptitle('Angular Velocity')
subplot(2,3,1), hold on, grid on
ylabel('$\omega_x \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(V_measured(1, :, :))*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(1, :, :))*180/pi,'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(V_ideal(1, :, :))*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,2), hold on, grid on
ylabel('$\omega_y \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(V_measured(2, :, :))*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(2, :, :))*180/pi,'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(V_ideal(2, :, :))*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,3), hold on, grid on
ylabel('$\omega_z \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(V_measured(3, :, :))*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, squeeze(ukf_V(3, :, :))*180/pi,'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, squeeze(V_ideal(3, :, :)*180/pi),'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)
legend('Measured','Filtered','Ideal','interpreter','latex','fontsize',15)

subplot(2,3,4), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_measured(:, 1)*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, EulAng_ukf(:, 1)*180/pi,'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, EulAng_ideal(:, 1)*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,5), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_measured(:, 2)*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, EulAng_ukf(:, 2)*180/pi,'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, EulAng_ideal(:, 2)*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)

subplot(2,3,6), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_measured(:, 3)*180/pi,'color',[0.80,0.80,0.80],'linewidth',5)
plot(t, EulAng_ukf(:, 3)*180/pi,'linewidth',2,'color',[0.00,0.00,0.00],'linewidth',1)
plot(t, EulAng_ideal(:, 3)*180/pi,'--','linewidth',2,'color',[0.85,0.33,0.10],'linewidth',2)


set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])

%% UKF Perf
%Position
figure()
subplot(2,3,1), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('x [km]','interpreter','latex');
plot(t, squeeze(ukf_g(1, 4, :))-squeeze(g_ideal(1, 4, :)),'linewidth',2)
plot(t,3*std_x(4,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(4,:),'--','linewidth',2);

subplot(2,3,2), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('y [km]','interpreter','latex');
plot(t, squeeze(ukf_g(2, 4, :))-squeeze(g_ideal(2, 4, :)),'linewidth',2)
plot(t,3*std_x(5,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(5,:),'--','linewidth',2);

subplot(2,3,3), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('z [km]','interpreter','latex');
h1=plot(t, squeeze(ukf_g(3, 4, :))-squeeze(g_ideal(3, 4, :)),'linewidth',2)
h2=plot(t,3*std_x(6,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
h3=plot(t,-3*std_x(6,:),'--','linewidth',2);
legend([h2 h3],'$+3\sigma$ bound','$-3\sigma$ bound','interpreter','latex','Fontsize',12,'location','southeast')

subplot(2,3,4), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_z \,[km/s]$','interpreter','latex');
plot(t, squeeze(ukf_V(4, :, :))-squeeze(V_ideal(4, :, :)),'linewidth',2)
plot(t,3*std_x(10,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(10,:),'--','linewidth',2);

subplot(2,3,5), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_y \,[km/s]$','interpreter','latex');
plot(t, squeeze(ukf_V(5, :, :))-squeeze(V_ideal(5, :, :)),'linewidth',2)
plot(t,3*std_x(11,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(11,:),'--','linewidth',2);

subplot(2,3,6), hold on, grid on
xlabel('t [s]','interpreter','latex');
ylabel('$v_z \,[km/s]$','interpreter','latex');
plot(t, squeeze(ukf_V(6, :, :))-squeeze(V_ideal(6, :, :)),'linewidth',2)
plot(t,3*std_x(12,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(12,:),'--','linewidth',2);


set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])

%Angular Velocity
figure()
% suptitle('Angular Velocity')
subplot(2,3,1), hold on, grid on
ylabel('$\omega_x \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(ukf_V(1, :, :))-squeeze(V_ideal(1, :, :)),'linewidth',2)
plot(t,3*std_x(7,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(7,:),'--','linewidth',2);

subplot(2,3,2), hold on, grid on
ylabel('$\omega_y \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, squeeze(ukf_V(2, :, :))-squeeze(V_ideal(2, :, :)),'linewidth',2)
plot(t,3*std_x(8,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(8,:),'--','linewidth',2);


subplot(2,3,3), hold on, grid on
ylabel('$\omega_z \, [deg/s]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
h1=plot(t, squeeze(ukf_V(3, :, :))-squeeze(V_ideal(3, :, :)),'linewidth',2)
h2=plot(t,3*std_x(9,:),'--','linewidth',2,'color',[0.47,0.67,0.19]);
h3=plot(t,-3*std_x(9,:),'--','linewidth',2);
legend([h2 h3],'$+3\sigma$ bound','$-3\sigma$ bound','interpreter','latex','Fontsize',12,'location','southeast')

subplot(2,3,4), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_ukf(:, 1)*180/pi- EulAng_ideal(:, 1)*180/pi,'linewidth',2)
plot(t,3*std_x(1,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(1,:)*180/pi,'--','linewidth',2);

subplot(2,3,5), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_ukf(:, 2)*180/pi- EulAng_ideal(:, 2)*180/pi,'linewidth',2)
plot(t,3*std_x(2,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(2,:)*180/pi,'--','linewidth',2);

subplot(2,3,6), hold on, grid on
ylabel('$\Theta_1 \, [deg]$','interpreter','latex')
xlabel('t [s]','interpreter','latex');
plot(t, EulAng_ukf(:, 3)*180/pi- EulAng_ideal(:, 3)*180/pi,'linewidth',2)
plot(t,3*std_x(3,:)*180/pi,'--','linewidth',2,'color',[0.47,0.67,0.19]);
plot(t,-3*std_x(3,:)*180/pi,'--','linewidth',2);

set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','xlim'),'xlim',[t(1) t(end)])

%% Plot Orbit
S1 = [squeeze(g_ideal(1:3,4,:)); V_ideal(4:6,:)];
S2 = [squeeze(ukf_g(1:3,4,:)); ukf_V(4:6,:)];
S3 = [squeeze(g_measured(1:3,4,:)); V_measured(4:6,:)];

Extent = 5;
%[Asteroid radii]Axes extent.

AXES= 0;
%[-]Enable/Disable asteroid axes 1/0 (ON/OFF)

Projection= 0;
%[-]Enable/Disable projection plot 1/0 (ON/OFF)

PlotOrbit(AXES,Extent,Asteroid,Projection,S1(1:3,:),'Ideal Motion',S1(1:3,:),'Filtered Motion',...
S3(1:3,:),'Raw Motion');
%[-]Plots the orbit in three (3) dimensions.
