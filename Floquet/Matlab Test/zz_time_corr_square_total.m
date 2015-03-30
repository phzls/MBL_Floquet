% Construct ZZ time correlation for a state
function time_corr_square = zz_time_corr_square_total(V)

d = length(V(1,:));

left_square = 0;
right_square = 0;
left_right_square = 0;
left_ave_right_ave = 0;

for n = 1:d
    state = V(:,n);
    [l,r,l_r] = zz_time_corr_square(state);
    left_square = left_square + l*l;
    right_square = right_square + r*r;
    left_right_square = left_right_square + l_r*l_r;
    left_ave_right_ave = left_ave_right_ave + l*r;
end

left_square = left_square / d;
right_square = right_square / d;
left_right_square = left_right_square / d;
left_ave_right_ave = left_ave_right_ave / d;

time_corr_square = left_right_square - left_square * right_square - left_ave_right_ave * left_ave_right_ave;