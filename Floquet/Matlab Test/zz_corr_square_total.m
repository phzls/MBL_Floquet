% Construct ZZ correlation for a state
function [corr_square,left,right, left_right] = zz_corr_square_total(V)

d = length(V(1,:));

corr_square = [];
left = [];
right = [];
left_right = [];

for n = 1:d
    state = V(:,n);
    [corr,l,r,l_r] = zz_corr_square(state);
    corr_square = [corr_square, corr];
    left = [left,l];
    right = [right,r];
    left_right = [left_right,l_r];
end

