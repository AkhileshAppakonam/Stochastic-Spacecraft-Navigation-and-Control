function R = expmso3(Tet)
% The matrix exponential for so(3) --> SO(3) (a.k.a Rodrigues formula)
%
% Author: Morad Nazari
tet = norm(Tet);
R = eye(3) + sine(tet) * crossm(Tet) + cosine(tet) * crossm(Tet)^2;

function cos_ = cosine(x)
n_max = 10; % number of terms in Taylor expansion
if x > 0.3
    cos_ = (1 - cos(x)) / x^2;
else
    cos_ = 0;
    for n = 0:n_max
        cos_ = cos_ + (-1)^n * x^(2*n) / factorial(2*(n+1));
    end
end

function sin_ = sine(x)
n_max = 10; % number of terms in Taylor expansion
if x > 0.2
    sin_ = sin(x) / x;
else
    sin_ = 0;
    for n = 0:n_max
        sin_ = sin_ + (-1)^n * x^(2*n) / factorial(2*n+1);
    end
end