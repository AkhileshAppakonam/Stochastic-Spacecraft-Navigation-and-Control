function [Psi, g_, V_, tet_] = f_sys(dt, ti, g0_, V0_, w_k, Xbar)

% NOISE
g0_ = Xbar*expm(vedge(w_k(1:6)))*g0_;
V0_ = V0_ + w_k(7:12);

% INTEGRATOR
dt2 = dt/2;
dt6 = dt/6;

tt = ti + dt2;
tet_ = norm(logSO3(g0_(1:3, 1:3)));
[~, dgdt, dVdt] = dynamicsTSE3(ti, g0_, V0_);
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

[psi, dg0t, dV0t] = dynamicsTSE3(ti + dt, g0t, V0t);

V0new = V0_ + dt6 * (dVdt + dV0t + 2*dVm);

omega0new = V0new(1:3);
dTet0new = omega0new * dt6; % dt changed to dt6
R0new = R0 * expmso3(dTet0new);

r0new = g0_(1:3, 4) + dt6 * (dgdt(1:3, 4) + dg0t(1:3, 4) + 2*drm);
g0new = [R0new    r0new;
    zeros(1, 3)       1];

% OUTPUT
Psi = psi;
V0_ = V0new;
g0_ = g0new;
g_ = g0_;
V_ = V0_;
end


