function [g_, V_, tet_,X_] = rk4_dynamics_B(dt, ti, tf, X0_, V0_)
global lt 
dt2 = dt/2;
dt6 = dt/6;

R0t = zeros(3, 3);            dRm = zeros(3, 3);           R0new = zeros(3, 3);
X0t = zeros(3, 1);               dXm = zeros(3, 1);              X0new = zeros(3, 1);
X0t = zeros(4, 4);            dgm = zeros(4, 4);           g0new = zeros(4, 4);
omega0t = zeros(3, 1);           d_omega_m = zeros(3, 1);        omega0new = zeros(3, 1);
dTet0t = zeros(3, 1);            d_dTet_m = zeros(3, 1);         dTet0new = zeros(3, 1);
g_ = zeros(4, 4, lt);         V_ = zeros(6, lt);           X_ = zeros(6, lt);
u_ = zeros(6, lt);              tet_ = zeros(1, lt);

indx = 0; % current time index
fprintf('time index: \n');
for t = ti: dt: (tf)
    indx = indx + 1;    fprintf('            %g / %g \n', indx, round(tf/dt))
    tt = t + dt2;
    tet_(1, indx) = norm(X0_(1:3));
    [dXdt, dVdt] = dynamics_B_NEW(t, X0_, V0_);
    V0t = V0_ + dt2 * dVdt;
    X0t = X0_ + dt2 * dXdt;
    
    [dX0t, dV0t] = dynamics_B_NEW(tt, X0t, V0t);
    V0t = V0_ + dt2 * dV0t;
    X0t = X0_ + dt2 * dX0t;
    
    [dXmold, dVmold] = dynamics_B_NEW(tt, X0t, V0t);
    V0t = V0_ + dt * dVmold;
    X0t = X0_ + dt * dXmold;
    
    dVm = dV0t + dVmold;
    dXm = dX0t + dXmold;
    
    [dX0t, dV0t] = dynamics_B_NEW(t + dt, X0t, V0t);
    V0new = V0_ + dt6 * (dVdt + dV0t + 2*dVm);
    X0new = X0_ + dt6 * (dXdt + dX0t + 2*dXm);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial condition for t + dt       &         augmenting the tensors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X0_ = X0new;
    g0_ = expm(vedge(X0_)); 
    V0_ = V0new;
    g_(:, :, indx) = g0_;
    V_(:, indx) = V0_;
    X_(:, indx) = X0_;
end