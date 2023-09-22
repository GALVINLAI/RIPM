function result = normOfStruct(struct1)
    % Initialize the result as zero
    result = 0;
    
    % Loop over each field and compute the sum of squares, then accumulate
    fields = fieldnames(struct1);
    for i = 1:length(fields)
        fieldName = fields{i};
        result = result + sum(sum(struct1.(fieldName).^2));
    end
    
    % Compute the square root to get the norm
    result = sqrt(result);
end

