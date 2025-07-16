function tensor = TT_reconstruct(tt_cores)
%TT_reconstruct(tt,n) 
%   Algorithm that reconstructs a N-th order tensor from the tensor train
%   decomposition
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
%   tensor (N-dimensional double):  N-th order tensor reconstructed from tt. 
%------ Implement your code below ------
tensor = squeeze(tt_cores{1});

for i = 2:size(tt_cores, 2)-1
    tensor = tensor * hidden_mode_n_matricization(tt_cores{i}, 1);
    tensor = reshape(tensor, size(tensor,1)*size(tt_cores{i},2),[]);
end
end