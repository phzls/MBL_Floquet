function left_spin = markov_x_left_spin(L,T,J,K, jump, Ur, state, debug)

tau = 0.8;
[U_up, H] = XXZ_bath_x_binary(L, J, K, tau, 1);
[U_down, H] = XXZ_bath_x_binary(L, J, K, tau, -1);

U_up = Ur * U_up;
U_down = Ur * U_down;

U_up

U_down

U_up_used = (U_up)^jump;
U_down_used = (U_down)^jump;

U_up_used

U_down_used


% Given initial state, compute density matrix
density = state * state';

if debug
    disp(0);
    disp(density);
end

% Left spin
left_spin = [];

% Time evolution
for t=0:T-1
    d = diag(density);
    left = 0;
    for i=1:length(d)
        spin = mod( idivide(int32(i-1), int32(2^(L-1))), 2 );
        spin = double(spin);
        if spin == 0
            left = left - d(i);
        else
            left = left + d(i);
        end
    end
    left_spin = [left_spin, left];
    temp = U_up_used * density * U_up_used' + U_down_used * density * U_down_used';
    density = 0.5*temp;
    if debug
        disp(t+1);
        disp(density);
    end
end


