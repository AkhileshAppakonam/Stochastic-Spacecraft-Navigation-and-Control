function [psi,gdot, Vdot] = dynamicsTSE3(t, g, V)

%Variables
global II control UU
R = g(1:3,1:3);
r = g(1:3,4);

%gdot
gdot = g * vedge(V);

%Pan Model
[F,M] = Pan_GravityModel(r,R,V,t);

%Control Input
if control=="yes"
[psi,u] = ControlSystem(g,V);
else 
u = zeros(6,1);  
psi = 0;  
end


%Vdot
Vdot = II \ ad_conj(V) * II * V  + II \ [M;F] + II \ (u);

UU = [UU u];
end