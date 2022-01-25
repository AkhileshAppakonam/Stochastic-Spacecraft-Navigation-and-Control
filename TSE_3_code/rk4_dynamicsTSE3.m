function [Psi, g_, V_, tet_] = rk4_dynamicsTSE3(dt, ti, tf, g0_, V0_)
global lt 
dt2 = dt/2;
dt6 = dt/6;

R0t = zeros(3, 3);            dRm = zeros(3, 3);           R0new = zeros(3, 3);
r0t = zeros(3, 1);               drm = zeros(3, 1);              r0new = zeros(3, 1);
g0t = zeros(4, 4);            dgm = zeros(4, 4);           g0new = zeros(4, 4);
omega0t = zeros(3, 1);           d_omega_m = zeros(3, 1);        omega0new = zeros(3, 1);
dTet0t = zeros(3, 1);            d_dTet_m = zeros(3, 1);         dTet0new = zeros(3, 1);
Psi = zeros(6, lt);
g_ = zeros(4, 4, lt);         V_ = zeros(6, lt);           X_ = zeros(6, lt);
u_ = zeros(6, lt);            tet_ = zeros(1, lt);         Z = zeros(12, 1, lt); 
indx = 0; % current time index
fprintf('time index: \n');
for t = ti: dt: (tf)
    
    indx = indx + 1;    fprintf('            %g / %g \n', indx, round(tf/dt))
    tt = t + dt2;
    tet_(1, indx) = norm(logSO3(g0_(1:3, 1:3)));
    [~, dgdt, dVdt] = dynamicsTSE3(t, g0_, V0_);
    V0t = V0_ + dt2 * dVdt;
    omega0t = V0t(1:3);
    dTet0t = omega0t * dt2; % dt changed to dt2
    R0 = g0_(1:3, 1:3, :);   
    R0t = R0 * expmso3(dTet0t);
    r0t = g0_(1:3, 4) + dt2 * dgdt(1:3, 4);
    g0t = [R0t    r0t;
        zeros(1,3)      1];

    [~, dg0t, dV0t] = dynamicsTSE3(tt, g0t, V0t);
    
    V0t = V0_ + dt2 * dV0t;
    
    omega0t = V0t(1:3);
    dTet0t = omega0t * dt2; % dt changed to dt2
    R0t = R0 * expmso3(dTet0t);
    %%%%%%%%%%%%%%%%%%%%%%%%
    r0t = g0_(1:3, 4) + dt2 * dg0t(1:3, 4);
    g0t = [R0t    r0t;
        zeros(1,3)      1];
    
    [~, dgmold, dVmold] = dynamicsTSE3(tt, g0t, V0t);
    V0t = V0_ + dt * dVmold;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    omega0t = V0t(1:3);
    dTet0t = omega0t * dt;
    R0t = R0 * expmso3(dTet0t);
    %%%%%%%%%%%%%%%%%%%%%%%%
    r0t = g0_(1:3, 4) + dt * dgmold(1:3, 4);
    g0t = [R0t    r0t;
        zeros(1,3)      1];
    
    dVm = dV0t + dVmold;
    
    d_omega_m = dVm(1:3);
    d_dTet_m = d_omega_m * dt;
    dRm = expmso3(d_dTet_m);
    
    drm = dg0t(1:3, 4) + dgmold(1:3, 4);
    dgm = [dRm            drm;
           zeros(1,3)      1];
    
    [psi, dg0t, dV0t] = dynamicsTSE3(t + dt, g0t, V0t);
    
    V0new = V0_ + dt6 * (dVdt + dV0t + 2*dVm);
    
    omega0new = V0new(1:3);
    dTet0new = omega0new * dt6; % dt changed to dt6
    R0new = R0 * expmso3(dTet0new);
    
    r0new = g0_(1:3, 4) + dt6 * (dgdt(1:3, 4) + dg0t(1:3, 4) + 2*drm);
    g0new = [R0new    r0new;
        zeros(1, 3)       1];    
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial condition for t + dt       &         augmenting the tensors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Psi(:, indx) = psi;
    V0_ = V0new;
    g0_ = g0new;
    g_(:, :, indx) = g0_;
    V_(:, indx) = V0_;
end




end


