function ent = markov(L,T,J,Ur, state, debug)

tau = 0.8;
[U_up, H] = XXZ_binary(L, J*tau, 1);
[U_down, H] = XXZ_binary(L, J*tau, -1);

U_up = Ur * U_up;
U_down = Ur * U_down;

% Given initial state, compute density matrix
density = state * state';

if debug
    disp(0);
    disp(density);
end

% Entropy
ent = [];

% Time evolution
for t=0:T-1
    V = eig(density);
    V = V.*log2(V);
    ent = [ent, -sum(V)];
    temp = U_up * density * U_up' + U_down * density * U_down';
    density = 0.5*temp;
    if debug
        disp(t+1);
        disp(density);
    end
end


