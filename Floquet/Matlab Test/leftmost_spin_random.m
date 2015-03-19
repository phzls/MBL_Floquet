% Construct a density matrix of a random state which has a random pure spin
% state at the leftmost and an identity density matrix for the rest
function rho = leftmost_spin_random(L, theta, phi)

% Identity matrix on each site
I = [1,0;0,1];

rho = I;
for n=2:L-1
    rho = kron(I,rho);
end

rho = rho / (2^(L-1));

leftmost_rho = [(1-cos(theta))/2, exp(j*phi)*sin(theta)/2; 
    exp(-j*phi)*sin(theta)/2, (1+cos(theta))/2 ];

rho = kron(leftmost_rho, rho);