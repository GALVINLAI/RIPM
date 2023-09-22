function minValue = findMinInStruct(struct1)
    % Initialize the minValue as positive infinity for comparison
    minValue = Inf;
    
    % Loop over each field and find the minimum value
    fields = fieldnames(struct1);
    for i = 1:length(fields)
        fieldName = fields{i};
        currentMin = min(struct1.(fieldName)(:)); % min of all elements in the current field
        if currentMin < minValue
            minValue = currentMin;
        end
    end
end

