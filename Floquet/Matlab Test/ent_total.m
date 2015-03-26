% Construct entropy for eigenstates
function ent = ent_total(V)

d = length(V(1,:));

ent = [];

for n = 1:d
    state = V(:,n);
    s = entropy(state);
    ent = [ent, s];
end

