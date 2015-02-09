% Construct the isolated xxz floquet in basic binary basis
function [U, H] = XXZ_iso_binary(L, tau)

h = 0.8090;
g = 0.9045;

% Identity matrix on each site
I = [1,0;0,1];
sigma_x = [0,1;1,0];
sigma_z = [-1,0;0,1];
rank = 2^L;

H = zeros(rank);

% Single site part for Hamiltonian. Positions should be counted from right
for n=1:L
    H = H + g * tensor_single(n,sigma_x,L);
    H = H + h * tensor_single(n,sigma_z,L);
end

% Nearest Neighbor Interaction
for n=1:L-1
    H = H + tensor_double(n,sigma_z,L);
end

% Construct time evolution
U = expm(-j*tau*H);
