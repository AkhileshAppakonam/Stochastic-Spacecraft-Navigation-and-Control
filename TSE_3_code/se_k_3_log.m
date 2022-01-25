function [xi] = se_k_3_log(chi)
%SE_K_3_LOG logarithm
%
% Syntax:  [xi] = se_k_3_log(chi)
%
% Inputs:
%    chi - matrix
%
% Outputs:
%    phi - vector

phi = so3_log(chi(1:3, 1:3));
Xi = so3_inv_left_jacobian(phi)*chi(1:3, 4:end);
xi = [phi;
     Xi(:)];
end

function [phi] = so3_log(Rot)
%SO3_LOG logarithm
%
% Syntax:  [phi] = so3_log(Rot)
%
% Inputs:
%    Rot - rotation matrix
%
% Outputs:
%    phi - vector

TOL = 1e-9;
cos_angle = 0.5 * trace(Rot) - 0.5;
% Clip cos(angle) to its proper domain to avoid NaNs from rounding errors
cos_angle = min(max(cos_angle, -1), 1);
angle = acos(cos_angle);

% If angle is close to zero, use first-order Taylor expansion
if norm(angle) < TOL
    phi = so3_vee(Rot - eye(3));
else
    % Otherwise take the matrix logarithm and return the rotation vector
    phi = so3_vee((0.5 * angle / sin(angle)) * (Rot - Rot'));
end
end

function [phi] = so3_vee(Phi)
%SO3_VEE vee operator
%
% Syntax:  [phi] = so3_vee(Phi)
%
% Inputs:
%    phi - vector
%
% Outputs:
%    Phi - matrix

phi = [Phi(3, 2);
       Phi(1, 3);
       Phi(2, 1)];
end

function [J] = so3_inv_left_jacobian(phi)
%SO3_INV_LEFT_JACOBIAN inverse of Jacobian
%
% Syntax:  [J] = so3_inv_left_jacobian(phi)
%
% Inputs:
%    phi - vector
%
% Outputs:
%    J - Jacobian

TOL = 1e-9;

angle = norm(phi);
if angle < TOL
    % Near |phi|==0, use first order Taylor expansion
    J = eye(3) - 1/2 * so3_wedge(phi);
else
    axis = phi / angle;
    half_angle = angle/2;
    J = half_angle * cot(half_angle) * eye(3) + ...
        (1 - half_angle * cot(half_angle)) * (axis*axis') -...
        half_angle * so3_wedge(axis);
end
end



function [Phi] = so3_wedge(phi)
%SO3_WEDGE Wedge operator
%
% Syntax:  [Phi] = so3_wedge(phi)
%
% Inputs:
%    phi - vector
%
% Outputs:
%    Phi - matrix

Phi = [0, -phi(3), phi(2);
       phi(3), 0, -phi(1);
       -phi(2), phi(1), 0];
end
