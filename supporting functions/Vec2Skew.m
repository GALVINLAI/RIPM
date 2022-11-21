function A = Vec2Skew(v,isometric)
% inverse of function: v = Skew2Vec(A)
% Converts a vector "v" to a skew-symmetric matrix "A".
% Use isometric isomorphism (keep inner product) , if isometric=true.

if ~exist('isometric', 'var')
    isometric = 0;
end

if isometric
    v=v/sqrt(2);
end

% compute size of skew-symmetric matrix
l = length(v);
n = (1 + sqrt(1 + 8*l))/2; 

% Focus on the bottom-left element, starting from the first column on the
% left and working from top to bottom in each column.
A = zeros(n); 
I = [0 cumsum(n-1:-1:1)];
for m = 1:n-1
      A(m+1:n, m) = v(I(m)+1:I(m + 1));
end

% Insert skew-symmetric part (Top right part)
for m = 1:n-1
      A(m, m+1:n) = -A(m+1:n, m).';
end

end