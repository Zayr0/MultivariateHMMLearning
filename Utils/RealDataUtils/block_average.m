function averagedData = block_average(data, numBlocks)
    averagedData = zeros(size(data, 1), numBlocks);
    
    for i = 1:numBlocks
        averagedData(:, i) = mean(data(:, numBlocks*(i-1)+1:numBlocks*i), 2);
    end
end