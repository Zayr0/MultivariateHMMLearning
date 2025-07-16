function [C,U]=mlsvd_Nd(T)

% MLSVD_ND computes the MLSVD of an input tensor of any order.
%
% INPUT:
%   T (N-D array): original N-Dimensional tensor 
%
% OUTPUT (in the respective order): 
%   C (N-D array): MLSVD core tensor
%   U (cell)     : MLSVD factor matrices of each mode stored in a cell,
%                  where U{1,n} gives the nth mode factor matrix
% 
% Remarks: 
%   You need to return the following variables:
N = ndims(T);
U = cell(1,N); % U{1,n} gives the nth mode factor matrix
C = T;


% % YOUR CODE GOES HERE
for n=1:N
    % Unfold tensor over dimension n
    Tn = mode_n_matricization(T, n);
    % Calculate svd of unfolding
    [u,~,~] = svd(Tn,"econ");
    U{n} = u;
    
    % Update iteration of the core tensor
    C = mode_n_product(C, u', n);
end
end