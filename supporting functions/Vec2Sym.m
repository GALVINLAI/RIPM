function A = Vec2Sym(V)

% Converts vector "V" to an N-by-N matrix "A" 
l = length(V); 
n = (-1 + sqrt(1 + 8*l))/2; 
A = zeros(n); 
I = [0 cumsum(1:n)];
for m = 1:n
      A(1:m, m) = V(I(m)+1:I(m + 1));
end
% Insert symmetric part
for m = 1:n-1
      A(m+1:n, m) = A(m, m+1:n).';

end

end