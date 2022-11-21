function B=rand_Nonneg_St(n,k)
% generate a random point x such that
% x>=0 and x'*x=I.
if n<k
    error('rand_Nonneg_St: we should have that n >= k.')
end

t=floor(n/k);
y=n-t*k;

M = spherefactory(t);
B=[];
for i=1:k-1
    B = blkdiag(B,abs(M.rand()));
end

M = spherefactory(t+y);
B = blkdiag(B,abs(M.rand()));
B = B(randperm(size(B, 1)), :);
end

