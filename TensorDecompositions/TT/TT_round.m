function tt = TT_round(tt,epsilon)
%TT_round (tt,n) 
%   Algorithm that rounds a TT algorithm up to a defined approximation
%   error
%      
%INPUT:
%   tt (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%   epsilon (double):               error-bound for approximation error of TT decomposition
%OUTPUT:
%   tt_cores (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%------ Implement your code below ------
tt = false; %Remove this, when you implement your code.
end