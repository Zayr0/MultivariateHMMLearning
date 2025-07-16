function outer = outer_product(varargin)
    % OUTER_PRODUCT takes a list of variable number of vectors varargin as input and
    % computes the outer product between them explicitly.
    % INPUT list of variable number of vectors varargin.
    % OUTPUT outer product tensor.
    
    outer = varargin{1};
    order = size(outer,1);
    
    for i = 2:length(varargin)
        vector = varargin{i};
        order(i)=size(vector,1);
        outer = outer*(vector');
        outer = reshape(outer,[],1);
    end

    outer = reshape(outer,order);
end
