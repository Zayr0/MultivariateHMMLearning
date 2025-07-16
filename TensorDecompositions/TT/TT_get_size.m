function sz = TT_get_size(tt)
%TT_get_size (tt,n) 
%   Algorithm that returns the size of a tt
%      
%INPUT:
%   tt (cell array with N+1 cells):  tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%OUTPUT:
%   sz (N x 3 double):              size of TT decomposition, see example
%                                   figure 
%------ Implement your code below ------

iter = numel(tt)-1;
sz = zeros(iter, 3);

for i = 1:iter

    if (ndims(tt{i}) == 3)
        sz(i, :) = size(tt{i});
    elseif ((ndims(tt{i})==2) && (size(tt{i}, 1)==1))
        sz(i, :) = [1, size(tt{i})];
    elseif ((ndims(tt{i})==2) && (size(tt{i}, 3)==1))
        sz(i, :) = [size(tt{i}), 1];
    else
        display("Error: incorrect size of entry of TT");
    end
end
end