% Construct density states for isolated inter_random_floquet with largest
% leftmost spin z
function [state, density, left_spin] = Largest_Leftmost_Spin_Z(L, J, tau, Ur)

% XXZ part
[Uxxz, Hxxz] = XXZ_iso_binary(L,J,tau);

U = Ur * Uxxz;


% Diagonalize the matrix in binary basis
[V,D] = eig(U)

% Record average leftmost spin z value
leftmost = [];

[row,col] = size(V);

for i=1:col
    left = 0;
    for j=1:row
        spin = mod( idivide(int32(j-1), int32(2^(L-1))), 2 );
        spin = double(spin);
        left = left + (2*spin-1)*(norm(V(j,i)))^2;
    end
    leftmost = [leftmost,left];
end

leftmost

% Sort according to leftmost
[leftmost,I] = sort(leftmost);

state = V(:,I(end));

% Given initial state, compute density matrix
density = state * state';

left_spin = leftmost(end);
        