function noe = TT_numel(tt)
%TT_numel(tt_1,tt_2,tt_sum) 
%  computes numbrt of elements of TT-decomposition
%      
%INPUT:
%   tt (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, which has been
%                                   computed with hidden_TT_SVD right
%                                   before. Here, the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
% %OUTPUT:
%   noe (double):                   number of elements of TT-decompositions, 
%                                   which are stored in first N cores. 
%                                   The element stored in the N+1th core
%                                   does NOT(!) count towards the number of
%                                   elements of TT-decomposition. 
%------ Implement your code below ------
noe = 0;

for i = 1:size(tt, 2)-1
    noe = noe + numel(tt{i});
end
end
