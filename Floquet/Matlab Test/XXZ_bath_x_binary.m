% Construct the xxz floquet with bath at the right in basic binary basis.
% The bath is coupled by sigma_x
function [U, B] = XXZ_bath_x_binary(L, J, K, tau, bath)

if bath ~= 1 && bath ~= -1
    disp(bath);
    error('bath must be either 1 or -1')
end

h = 0.8090;
g = 0.9045;

% Identity matrix on each site
I = [1,0;0,1];
sigma_x = [0,1;1,0];
sigma_z = [-1,0;0,1];
rank = 2^L;

Hx = zeros(rank);
Hz = zeros(rank);

% Single site part for Hamiltonian. Positions should be counted from right
for n=1:L
    Hx = Hx + J*g * tensor_single(n,sigma_x,L);
    Hz = Hz + J*h * tensor_single(n,sigma_z,L);
end

% Bath part
B = bath*tensor_single(1, sigma_x, L);

% Nearest Neighbor Interaction
for n=1:L-1
    Hz = Hz + J*tensor_double(n,sigma_z,L);
end

B = expm(-j*K*tau*B);

% Construct time evolution
U = expm(-j*tau*Hx)*expm(-j*tau*Hz)*B;
