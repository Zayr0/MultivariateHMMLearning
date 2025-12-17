function data = loadRealData(DataEnvelope, numSamplesPerBlock, numBlocks, numChannels)
    data = zeros(3, numSamplesPerBlock * numBlocks);
    
    for i = 1:numChannels
        for j = 1:numBlocks-1
            data(i, ((j-1)*numSamplesPerBlock+1):numSamplesPerBlock*j) = DataEnvelope.(i){j};
        end
    end
end