% Given a state, compute its leftmost spin value
function left_spin = Leftmost_spin_z(state)
[row,col] = size(state);
L = log2(row);
left_spin = 0;
for j=1:row
    spin = mod( idivide(int32(j-1), int32(2^(L-1))), 2 );
    spin = double(spin);
    left_spin = left_spin + (2*spin-1)*(norm(state(j)))^2;
end