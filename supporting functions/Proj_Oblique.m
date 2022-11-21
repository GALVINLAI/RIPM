function Z=Proj_Oblique(X)
% Projection onto oblique manifold whose size is euqal to X.
% The oblique manifold OB(n,m) (the set of matrices of size nxm with
% unit-norm columns)

[n,k]=size(X);
Z=zeros(n,k);
for j=1:k
Z(:,j)=X(:,j)/norm(X(:,j));
end
end