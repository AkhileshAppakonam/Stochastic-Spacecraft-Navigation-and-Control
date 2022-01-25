function log_SO3 = logSO3(C)
% For the matrix logarithm on SO(3) we use the
%    "Proportional Derivative (PD) control on the Euclidean Group"
%    by F. Bullo and R.M. Murray.
% OR
%    expression (5.9) in "Geometric Control of Mechanical Systems"
%    by F. Bullo and A. D. Lewis.
%
% Modifications made by Amit K. Sanyal in Aug 2011 to account for
%   rotations by an odd multiple of pi radians.
% 
% More modification and comments made by Morad Nazari as well as 
%   constructing the Bernoulli numbers and Taylor expansion of x*csc(x) (April 2014)

tet = acos((trace(C)-1)/2);

nA = norm(C-C', 'fro'); % Frobenius norm of C-C'
E = (C - eye(3)) / 2;
ed = diag(E);
nE = norm(E, 'fro');
[med, i] = max(ed);

if nA < 1e-5 && nE < 1e-5
    log_SO3 = zeros(3,1);  % C ~ C' (small rotations --> rotation matrix is almost identity)
elseif nA < 1e-5 && nE > 0.1
    ri = sqrt(1+med);
    if i == 1
        log_SO3 = pi * [ri; E(1,2)/ri; E(1,3)/ri];
    elseif i == 2
        log_SO3 = pi * [E(2,1)/ri; ri; E(2,3)/ri];
    else log_SO3 = pi * [E(3,1)/ri; E(3,2)/ri; ri];
    end
else
    logC = xcsc(tet) * (C-C');
    log_SO3 = unskew(logC);
end

function xcsc_ = xcsc(x)
% constructing the x*csc(x) function
n_max = 20;
if abs(x) < 0.2
    % constructing the Bernoulli numbers (note that B0 = 1)
    B = zeros(1, 2*n_max);
    for n = 1: 2*n_max
        B(n) = - 1/(n+1); % first term in the Bernoulli series accounting for m = 0
        for m = 1: n-1
            B(n) = B(n) - nchoosek(n,m) * B(m) / (n-m+1);
        end
    end
    xcsc_ = 1; % 1 is the first term of the Taylor expantion of x*csc(x)
    for n = 1: n_max
        xcsc_ = xcsc_ + B(2*n) * (-1)^n * (2-4^n) * x^(2*n) / factorial(2*n); 
    end
    xcsc_ = xcsc_ / 2; %  Taylor expantion of x*csc(x)/2
else
    xcsc_ = x * csc(x) / 2;
end
