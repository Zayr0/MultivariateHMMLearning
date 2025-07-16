function outer = outer_product_kron(varargin)
    % OUTER_PRODUCT_KRON takes a list of variable number of vectors varargin as input and
    % computes the outer product between them using the kron function which
    % is builtin in MATLAB.
    % INPUT list of variable number of vectors varargin.
    % OUTPUT outer product tensor.
    
    outer = varargin{end};
    order(length(varargin)) = size(outer,1);

    for i=(length(varargin)-1):-1:1 
        outer = kron(outer,varargin{i});

        order(i) = size(varargin{i},1);
    end
    outer = reshape(outer,order);
end
