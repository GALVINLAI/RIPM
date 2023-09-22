function result = multiplyStructWithScalar(struct1, scalar)
    % Initialize the result structure
    result = struct();
    
    % Loop over each field and multiply every element by the scalar
    fields = fieldnames(struct1);
    for i = 1:length(fields)
        fieldName = fields{i};
        result.(fieldName) = struct1.(fieldName) * scalar;
    end
end

