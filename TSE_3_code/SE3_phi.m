function [new_state] = SE3_phi(state, xi)
% Inputs:
%    state - state
%    xi - uncertainty
%
% Outputs:
%    new_state - state

new_state.g = state.g * expm(vedge(xi(1:6)));
new_state.V = state.V + xi(7:12);
end