% Construct the xxz floquet with bath at the right in basic binary basis.
% The bath is coupled by sigma_x
function [U, B] = XXZ_random_bath_x_binary(L, lambda, K, random_h, tau, bath)

if bath ~= 1 && bath ~= -1
    disp(bath);
    error('bath must be either 1 or -1')
end

[U_iso, Hx, Hz] = XXZ_random_binary(L, lambda, random_h, tau);

% Identity matrix on each site
I = [1,0;0,1];
sigma_x = [0,1;1,0];
sigma_z = [-1,0;0,1];
rank = 2^L;

% Bath part
B = bath*tensor_single(1, sigma_x, L);

B = expm(-j*K*tau*B);

% Construct time evolution
U = U_iso*B;
