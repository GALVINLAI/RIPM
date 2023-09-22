function result = maxZeroStruct(struct1)
    % Initialize the result structure
    result = struct();
    
    % Loop over each field and apply the max(,0) operation element-wise
    fields = fieldnames(struct1);
    for i = 1:length(fields)
        fieldName = fields{i};
        result.(fieldName) = max(struct1.(fieldName), 0);
    end
end

