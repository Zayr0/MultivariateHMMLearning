function ranks = TT_get_ranks(tt)
%TT_get_ranks (tt,n) 
%   Algorithm that returns the ranks of a tt
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
%   ranks (N+1 x 1 double):         ranks of TT decomposition, see example
%                                   figure 
%------ Implement your code below ------

iter = length(tt)-1;
ranks = zeros(iter, 1);
for i= 1:iter
    if (i == iter)
        ranks(i) = size(tt{i}, 2);
    else
        ranks(i) = size(tt{i}, 3);
    end
end
end