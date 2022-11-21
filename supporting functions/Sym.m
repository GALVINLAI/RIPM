function result=Sym(A)
    % projector onto symmetric matrix space
    result=0.5*(A+A.');
end