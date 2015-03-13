% Construct the isolated random xxz floquet in basic binary basis
function [U, Hx, Hz] = XXZ_random_binary(L, lambda, random_h, tau)

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
    Hx = Hx + (1-lambda)*g * tensor_single(n,sigma_x,L);
    Hz = Hz + (h + lambda*random_h(n))* tensor_single(n,sigma_z,L);
end

% Nearest Neighbor Interaction
for n=1:L-1
    Hz = Hz + tensor_double(n,sigma_z,L);
end

% Construct time evolution
U = expm(-j*tau*Hx)*expm(-j*tau*Hz);
