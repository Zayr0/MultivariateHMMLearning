function tt = site_n(tt,n)
%site_n(tt,n) 
%   Algorithm that brings tensor train into site_n mixed canoncical form. 
%   decomposition 
%      
%INPUT:
%   tt (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0. 
%   n (int):                        n with 1<=n<=N, location where norm is 
%                                   supposed to be moved to.
%OUTPUT:
%   tt (cell array with N+1 cells): tensor train (tt) decomposition with 
%                                   norm in positin n. It is a tt of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell.
%------ Implement your code below ------
tt = false; %Remove this, when you implement your code.
end