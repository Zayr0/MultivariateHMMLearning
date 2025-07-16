function Sum  = TT_add(A,B)
%TT_add(tt,n) 
%   Algorithm adds two tensors in tensor train (TT) format
%      
%INPUT:
%   A (cell array with N+1 cells):  tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%   B (cell array with N+1 cells):  tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%OUTPUT:
%   Sum (cell array with N+1 cells): sum of the TTs A and B in TT format. 
%                                   The TT-cores are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell.If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%------ Implement your code below ------
Sum = false;%Remove this, when you implement your code.
end