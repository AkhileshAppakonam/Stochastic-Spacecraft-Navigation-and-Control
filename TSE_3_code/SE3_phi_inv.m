function [xi] = SE3_phi_inv(state, hat_state)
% Inputs:
%    state - state
%    hat_state - state
%
% Outputs:
%    xi - uncertainty

chi = [state.g(1:3,1:3) state.g(1:3,4);
       zeros(1, 3), 1];
hat_chi = [hat_state.g(1:3,1:3) hat_state.g(1:3,4);
           zeros(1, 3), 1];

xi = [se_k_3_log((se_k_3_inv(hat_chi) * chi));
      state.V- hat_state.V];

% chi = [state.g(1:3,1:3) state.g(1:3,4) state.V(1:3), state.V(4:6);
%        zeros(3, 3), eye(3,3)];
% hat_chi = [hat_state.g(1:3,1:3) hat_state.g(1:3,4)  hat_state.V(1:3), hat_state.V(4:6);
%        zeros(3, 3), eye(3,3)];
% 
% xi = se_k_3_log((se_k_3_inv(hat_chi) * chi));


end

function [chi_inv] = se_k_3_inv(chi)
% Inputs:
%    chi - matrix
%
% Outputs:
%    chi_inv - matrix

k = length(chi) - 3;
chi_inv = [chi(1:3, 1:3)' -chi(1:3, 1:3)'*chi(1:3, 4:end);
    zeros(k, 3) eye(k)];
end