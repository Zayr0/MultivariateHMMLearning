function B_matched = match_states(A, B)
    % Ensure A and B are column vectors
    A = A(:);
    B = B(:);
    
    unique_B = unique(B);
    B_mapped = zeros(size(B));
    
    % For each unique value in B, find what it most often corresponds to in A
    mapping = containers.Map('KeyType', 'double', 'ValueType', 'double');

    for i = 1:length(unique_B)
        val = unique_B(i);
        indices = (B == val);
        corresponding_A = A(indices);

        if isempty(corresponding_A)
            continue;
        end

        % Majority value in A corresponding to this B state
        mode_val = mode(corresponding_A);
        mapping(val) = mode_val;
    end

    % Apply the mapping
    for i = 1:length(B)
        B_mapped(i) = mapping(B(i));
    end

    B_matched = B_mapped;
end
