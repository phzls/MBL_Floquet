% Compute half-chain entanglement entropy given a state
function ent = entropy(state)
L = int32(log2(length(state)));
left = L/2;

density = zeros(2^(left));
for n = 0:2^(left)-1
    for m1 = 0:2^(left)-1
        for m2 = 0:2^(left)-1 
            num1 = n + m1 * 2^(left);
            num2 = n + m2 * 2^(left);
            density(m1+1,m2+1) = density(m1+1,m2+1) + state(num1+1) * conj(state(num2+1));
        end
    end
end

D = eig(density);
ent = 0;
for n=1:length(D)
    if abs(D(n)) > 1.0e-15
        ent = ent - D(n) * log2(D(n));
    end
end

            
       

