% Compute time evolution of a leftmost spin random state under xxz random
% floquet operator
function [spin_x, spin_y, spin_z] = leftmost_spin_xxz_random_evolution(L, lambda, random_h, tau, theta, phi, T, debug)

tau = 0.8

[U, Hx, Hz] = XXZ_random_binary(L, lambda, random_h, tau)

% Identity matrix on each site
I = [1,0;0,1];
sigma_x = [0,1;1,0];
sigma_z = [-1,0;0,1];
sigma_y = [0,j;-j,0];

id = I
for n=2:L-1
    id = kron(I,id);
end

left_x = kron(sigma_x, id);
left_y = kron(sigma_y, id);
left_z = kron(sigma_z, id);

density = leftmost_spin_random(L, theta, phi);

if debug
    disp(0);
    disp(density);
end

% Left spin
spin_x = [];
spin_y = [];
spin_z = [];

% Time evolution
for t=0:T-1
    x_value = trace(density*left_x);
    spin_x = [spin_x, x_value];
    y_value = trace(density*left_y);
    spin_y = [spin_y, y_value];
    z_value = trace(density*left_z);
    spin_z = [spin_z, z_value];
    
    density = U * density * U';

    if debug
        disp(t+1);
        disp(density);
    end
end