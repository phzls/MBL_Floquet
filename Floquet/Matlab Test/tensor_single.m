% Construct a matrix operating on the whole with a matrix appearing at
% position pos. The position is counted from right.
function A = tensor_single(pos, B, L)
if pos<1
    error('Position must be no less than 1')
end
I = [1,0;0,1];
if pos==1
    A = B;
else
    A = I;
    for i=2:pos-1
        A = kron(I,A);
    end
    A = kron(B,A);
end 

for i=pos+1:L
    A = kron(I,A);
end
end