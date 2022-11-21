function Z=Proj_Nonnega(X)
Z=X;
Z(Z<0)=0;
end