function Pbest = MultivariatePermutationFitEmission(Os, O_hats)
    
    n = length(Os);
    Pbests = cell(n, 1);
    
    for i = 1:n
        Pbests{i} = PermutationFitEmission(Os{i}, O_hats{i});
    end
    
    if(~isequal(Pbests{:}))
        warning("Not all permutations are the same")
        Pbest = findMostFrequentMatrix(Pbests);
    else
        Pbest = Pbests{1};
    end
end


function mostFrequentMatrix = findMostFrequentMatrix(cellArray)
    n = numel(cellArray);
    keys = strings(n, 1);
    map = containers.Map('KeyType', 'char', 'ValueType', 'any');

    for i = 1:n
        % Convert matrix to string to use as key
        key = mat2str(cellArray{i});
        keys(i) = key;

        if isKey(map, key)
            map(key) = [map(key), i];
        else
            map(key) = i;
        end
    end

    % Find the key with the most occurrences
    maxCount = 0;
    mostFrequentKey = '';
    for k = keys'
        key = char(k);
        if isKey(map, key)
            count = numel(map(key));
            if count > maxCount
                maxCount = count;
                mostFrequentKey = key;
            end
        end
    end

    % Return one instance of the most frequent matrix
    idx = map(mostFrequentKey);
    mostFrequentMatrix = cellArray{idx(1)};
end