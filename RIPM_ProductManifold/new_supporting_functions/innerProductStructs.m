function result = innerProductStructs(struct1, struct2)
    % Ensure the two structures have the same fields
    if ~isequal(fieldnames(struct1), fieldnames(struct2))
        error('The two structures do not have the same fields.');
    end
    
    % Initialize the result as zero
    result = 0;
    
    % Loop over each field and compute the inner product, then accumulate
    fields = fieldnames(struct1);
    for i = 1:length(fields)
        fieldName = fields{i};
        result = result + sum(sum(struct1.(fieldName) .* struct2.(fieldName)));
    end
end