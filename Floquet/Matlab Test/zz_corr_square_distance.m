% Construct ZZ correlation for a state for a given starting position of
% left sigma_z matrix
function [corr_square,left,right, left_right] = zz_corr_square_distance(state,d)
density = state * state';

% Identity matrix on each site
I = [1,0;0,1];
sigma_x = [0,1;1,0];
sigma_z = [-1,0;0,1];
sigma_y = [0,j;-j,0];


L = int32(log2(length(state)));

if d == 1
    right_z = sigma_z;
    left_z = sigma_z;
else
    right_z = I;
    left_z = I;
end

for n=2:d-1
    right_z = kron(I, right_z);
    left_z = kron(left_z,I);
end

if d>1
    right_z = kron(sigma_z, right_z);
    left_z = kron(left_z, sigma_z);
end

left_right_z = right_z;

for n=d+1:L-d
    left_right_z = kron(I, left_right_z);
end

left_right_z = kron(left_z, left_right_z);

for n=d+1:L
    right_z = kron(I, right_z);
    left_z = kron(left_z, I);
end


left = trace(density*left_z);
right = trace(density*right_z);
left_right = trace(density*left_right_z);

corr_square = (left_right - left*right)^2;