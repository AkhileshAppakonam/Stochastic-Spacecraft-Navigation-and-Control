clear var; clc; close all;
tic
dbstop if error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code simulates the control strategies given at the end of the slide-
% set named "Control Laws for a Single Rigid Body on SE(3)" plus a locally
% stabilizing control strategy proposed here.
%
% ------------------------------------------------------------------------
% R: rotation matrix from body frame to inertial frame
% r: position vector of the center of mass
% g = [R r; zeros(1,3) 1]:  configuration (pose) of the rigid body
% ------------------------------------------------------------------------
% definition of the abbreviations used for the controller:
% LeSaBu15JGCD:   D. Lee, A. Sanyal, E.A. Butcher, and D. Scheeres, JGCD, 2015
% LeSaBuSc14AST:  D. Lee, A. Sanyal, E.A. Butcher, and D. Scheeres, AST, 2014
% LeSaBuSc14TAES: D. Lee, A. Sanyal, E.A. Butcher, and D. Scheeres, IEEE, TAES, 2014
% Na16local:      Morad's local stabilizing control law
%
% All these control strategies use exponential coordinates
% The first three controllers above are almost globally asymptotically
% stabilizing while the last one is locally asymptotically stabilizing.
% ------------------------------------------------------------------------
% NOTES: 1) V in this code is the same as Omega in the slides mentioned above.
%        2) R = g(1:3, 1:3, :);  r = g(1:3, 4, :)
%        3) To reproduce the results as in the slides, use dt = 1 for the 
%           controllers that use exponential coordinates and use dt = 0.1
%           for the TSE(3) controller.
%        4) last run as of 10-28-16: dt =.1; tf=500; tracking only for
%        LeSaBu15JGCD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Morad Nazari and Juan Cuadra [October 2016]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global cont_choose pw lit tumbling II lt N K K1 K2 kappa Rz1 Ry Rz2 alpha beta gamma_ q e a a_ b_ c_
%% input
J = diag([4.97 6.16 8.37]);
%[kg.m^2]Inertia tensor

m = 60;
II = [J         zeros(3);
    zeros(3)  m*eye(3)];

% a_ = 20;   b_ = .02;    c_=.1; % plotting the desired path for the LeSaBu15JGCD control
a_ = 2;   b_ = .02;    c_=.1; % plotting the desired path for the LeSaBuSc14AST control

%% timespan and initial conditions
N = 1;

dt = .1; % OR dt = 1 (See the NOTE #3 in the comments above)
ti = 0; % initial time
tf = 400; % (sec) final time
t = ti: dt: tf;
lt = length(t);

psi = [pi/5; pi/5; pi/5];
Rz1 = [cos(psi(3)) -sin(psi(3)) 0;              sin(psi(3))  cos(psi(3)) 0;            0       0       1      ];
Ry =  [cos(psi(2))    0    sin(psi(2));             0             1      0;   -sin(psi(2))     0   cos(psi(2))];
Rz2 = [cos(psi(1)) -sin(psi(1)) 0;              sin(psi(1))  cos(psi(1)) 0;            0       0       1      ];
R0 = Rz2 * Ry * Rz1;
r0 = 10*[1;0;0];
g0 = [R0 r0; zeros(1,3) 1];
omega0 = [.2; -.3; .4]; v0 = [1; -.2; -.3];
V0 = [omega0; v0];
X0 = vedge_inv( logSE3( g0 ) ); % X here is like Xrel in tracking control


%% choice of the control law
cont_choose = questdlg('step1: Select the control from:', 'cont_choose', 'literature', ...
    'present work', 'literature');
switch cont_choose
    case 'present work'
        pw = questdlg('step2: Select the control from:', 'present work', 'TSE(3)', ...
            'exponential coordinates', 'TSE(3)');
        switch pw
            case 'exponential coordinates'
                % K = 2; alpha = 1; beta = .1;
                K = diag([2,2,2,1,1,1]); alpha = diag([1,1,1,.03,.03,.03]); beta = .1;
            case 'TSE(3)'                
                tumbling = questdlg(' ', 'tumbling', 'tumbling', ...
                                      'no tumbling', 'tumbling');
                switch tumbling
                    case 'tumbling'
                        psi = [pi/2; -pi/2; pi/2];
                        omega0 = 100*[.2; -.3; .4]; 
                    case 'no tumbling'
                        psi = [pi/5; pi/5; pi/5];
                        omega0 = [.2; -.3; .4]; 
                end
                kappa = .02;
                k11 = .2; k12 = .1; k21 = 1; k22 = .01;
                K1 = blkdiag(k11*eye(3), k12*eye(3));
                K2 = blkdiag(k21*eye(3), k22*eye(3));
                a = [1.2; 1.1; 1]; % a1 > a2 > a3 >= 1
                e = eye(3);
                
                Rz1 = [cos(psi(3)) -sin(psi(3)) 0;              sin(psi(3))  cos(psi(3)) 0;            0       0       1      ];
                Ry =  [cos(psi(2))    0    sin(psi(2));             0             1      0;   -sin(psi(2))     0   cos(psi(2))];
                Rz2 = [cos(psi(1)) -sin(psi(1)) 0;              sin(psi(1))  cos(psi(1)) 0;            0       0       1      ];
                R0 = Rz2 * Ry * Rz1;
                g0 = [R0 r0; zeros(1,3) 1];
                V0 = [omega0; v0];
        end
    case 'literature'
        lit = questdlg('Select the control law', 'literature', 'LeSaBu15JGCD', ...
            'LeSaBuSc14AST', 'LeSaBuSc14TAES', 'LeSaBuSc14AST');
        switch lit
            case 'LeSaBu15JGCD'
                alpha1 = .2; alpha2 = .6;
                alpha = blkdiag(alpha1*eye(3), alpha2*eye(3));
                beta1 = 2; beta2 = .7;
                beta = blkdiag(beta1*eye(3), beta2*eye(3));
            case 'LeSaBuSc14AST'
                % alpha = .083;
                % gamma_= .004;
                alpha = .03;
                gamma_= .1;
            case 'LeSaBuSc14TAES'
                alpha1 = .3; alpha2 = .1;
                alpha = blkdiag(alpha1*eye(3), alpha2*eye(3));
                gamma_ = .15;
                q = 1.5; % (1 < q < 2)
        end
end

%% integration
if strcmp(cont_choose, 'present work') == 1 && strcmp(pw, 'TSE(3)') == 1
    [psi, g, V, tet] = rk4_dynamicsTSE3(dt, ti, tf, g0, V0);
else
    [g, V, tet] = rk4_dynamics_B(dt, ti, tf, X0, V0);
end

%% plot
if strcmp(cont_choose, 'present work') == 1 && strcmp(pw, 'TSE(3)') == 1
    figure
    plot(t, squeeze(psi(1,:)), 'k', t, squeeze(psi(2,:)), '--r', t, squeeze(psi(3,:)),'.-b');
    xlabelentries = {'$$t$$ (s)'};
    xlabel(xlabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
    ylabelentries = {'$$\psi$$'};
    ylabel(ylabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
    grid;
end

figure
subplot(2,1,1);
plot(t, tet.*180/pi, 'k')
xlabelentries = {'$$t$$ (s)'};
xlabel(xlabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
ylabelentries = {'$$\theta$$ (deg)'};
ylabel(ylabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
grid;
subplot(2,1,2);
plot(t, squeeze(g(1, 4, :)), 'k', t, squeeze(g(2, 4, :)), 'r', t, squeeze(g(3, 4, :)), 'b', 'markersize', 14)
xlabelentries = {'$$t$$ (s)'};
xlabel(xlabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
ylabelentries = {'$$r$$ (m)'};
ylabel(ylabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
grid

% plots of velocities in the ECI frame
figure
subplot(2,1,1);
plot(t, squeeze(V(1, :)), 'k', t, squeeze(V(2, :)), '--r', t, squeeze(V(3, :)), '-.g', 'linewidth', 2);
xlabelentries = {'$$t$$ (sec)'};
xlabel(xlabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
ylabelentries = {'$$\omega$$ (rad/s)'};
ylabel(ylabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
grid on
subplot(2,1,2);
plot(t, squeeze(V(4, :)), 'k', t, squeeze(V(5, :)), '--r', t, squeeze(V(6, :)), '-.g', 'linewidth', 2);
xlabelentries = {'$$t$$ (sec)'};
xlabel(xlabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
ylabelentries = {'$$v$$ (m/s)'};
ylabel(ylabelentries, 'interpreter', 'latex', 'FontName', 'Times', 'fontsize', 14)
grid on

figure (10)
% scatter3(0,0,0)
hold on;
%x = a_*sin(20*b_*pi*t); y = a_*cos(20*b_*pi*t); z = 20*c_*t;
%u = gradient(x); v = gradient(y); w = gradient(z);  
%scale = 0;
%quiver3(x,y,z,u,v,w,scale); % plot3(rd = [sin(2*pi*t); cos(2*pi*t); t];)
%plot3(x,y,z); % plot3(rd = [sin(2*pi*t); cos(2*pi*t); t];)
plot3(0,0,0,'.');
line([0 9], [0 0], [0 0]);  
line([0 0], [0 9], [0 0]);
line([0 0], [0 0], [0 9]);

%%
tt=size(t,2);

BB = [g(1, 4, :) g(2, 4, :) g(3, 4, :)]; % Position data
A = g(1:3, 1:3, :) ;                     % Orientation data

myvizualization_mod(t,BB,A)

toc

