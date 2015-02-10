function left_spin = markov_left_spin(L,T,J,Ur, state, debug)

tau = 0.8;
[U_up, H] = XXZ_bath_binary(L, J, tau, 1);
[U_down, H] = XXZ_bath_binary(L, J, tau, -1);

U_up = Ur * U_up;
U_down = Ur * U_down;

U_up

U_down

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
    temp = U_up * density * U_up' + U_down * density * U_down';
    density = 0.5*temp;
    if debug
        disp(t+1);
        disp(density);
    end
end


