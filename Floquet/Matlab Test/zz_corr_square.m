% Construct ZZ correlation for a state
function [corr_square,left,right, left_right] = zz_corr_square(state)
density = state * state';

% Identity matrix on each site
I = [1,0;0,1];
sigma_x = [0,1;1,0];
sigma_z = [-1,0;0,1];
sigma_y = [0,j;-j,0];

id = I;
id_small = sigma_z;

L = int32(log2(length(state)));

for n=2:L-1
    id = kron(I,id);
end

for n = 2:L-1
    id_small = kron(I,id_small);
end

left_z = kron(sigma_z, id);
right_z = kron(id, sigma_z);
left_right_z = kron(sigma_z, id_small);

left = trace(density*left_z);
right = trace(density*right_z);
left_right = trace(density*left_right_z);

corr_square = (left_right - left*right)^2;