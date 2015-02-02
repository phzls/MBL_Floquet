% Construct a matrix operating on the whole with a matrix appearing at
% position pos and pos+1. Positions are counted from right.
function A = tensor_double(pos, B, L)
if pos<1
    error('Position must be no less than 1')
end
I = [1,0;0,1];
if pos>=L
    error('Position must be less than chain size')
end

if pos==1
    A = B;
else
    A = I;
    for i=2:pos-1
        A = kron(I,A);
    end
    A = kron(B,A);
end 

A = kron(B,A);

for i=pos+2:L
    A = kron(I,A);
end
end