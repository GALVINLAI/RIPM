function result = multiplyStructs(struct1, struct2)
    % Ensure the two structures have the same fields
    if ~isequal(fieldnames(struct1), fieldnames(struct2))
        error('The two structures do not have the same fields.');
    end
    
    % Initialize the result structure
    result = struct();
    
    % Loop over each field and multiply element-wise
    fields = fieldnames(struct1);
    for i = 1:length(fields)
        fieldName = fields{i};
        result.(fieldName) = struct1.(fieldName) .* struct2.(fieldName);
    end
end