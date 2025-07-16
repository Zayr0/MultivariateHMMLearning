function Z = khatri_rao(varargin)
    % KHATRI_RAO takes a list of matrices and returns the (right)
    % Khatri-Rao product.
    % INPUT list of variable number of matrices varargin.
    % OUTPUT outer product outer.
    
    Z = varargin{1};

    for i=2:length(varargin)
        temp = [];

        for j=1:size(Z,1)
            temp = [temp; Z(j,:).*varargin{i}];
        end

        Z = temp;
    end
end