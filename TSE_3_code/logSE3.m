function log_ = logSE3(R)
% input: relative configurations in SE(3) (current and desired)
% output: logarithm map, i.e. Xvedge in se(3)
%
% Author : Morad Nazari, April 2014

Tet = logSO3(R(1:3,1:3)); % principal rotation vector Theta
tet = norm(Tet); % principal rotation angle theta

S = eye(3,3) + cosine(tet)*crossm(Tet) + sine(tet)*crossm(Tet)^2;

log_ = [crossm(Tet)   S\R(1:3,4); 
            zeros(1,4)];

function cos_ = cosine(x)
n_max = 10; % number of terms in Taylor expansion
if x > 0.4
    cos_ = (1-cos(x))/x^2;
else
    cos_ = 0;
    for n = 0:n_max
        cos_ = cos_ + (-1)^n * x^(2*n) / factorial(2*(n+1));
    end
end

function sin_ = sine(x)
n_max = 10; % number of terms in Taylor expansion
if x > 0.3
    sin_ = (x-sin(x))/x^3;
else
    sin_ = 0;
    for n = 0:n_max
        sin_ = sin_ + (-1)^n * x^(2*n)/factorial(2*n+3);
    end
end