function [C,U1,U2,U3]=mlsvd_3d(T)


% MLSVD_3D computes MLSVD for 3rd order tensors.
%
% INPUT:
%   T (3D array): original 3-D tensor (of size I_1 x I_2 x I_3)
%
% OUTPUT (in the respective order): 
%   C (3D array): MLSVD core tensor (has the same size as T)
%   U1 (matrix) : MLSVD factor matrix of the first mode 
%                 (size: I_1 x I_1)
%   U2 (matrix) : MLSVD factor matrix of the second mode 
%                 (size: I_2 x I_2)
%   U3 (matrix) : MLSVD factor matrix of the third mode 
%                 (size: I_3 x I_3)
% 
% Remarks: 
%   Note that the computation of MLSVD is exact and deterministic. 


%% PART 1: Compute the factor matrices

% Hint: You can compute each n-mode factor matrix (U1, U2 and U3) 
% directly using the SVD of the n-mode unfolding of the tensor. 

% % YOUR CODE GOES HERE

T1 = mode_n_matricization(T, 1);
T2 = mode_n_matricization(T, 2);
T3 = mode_n_matricization(T, 3);


[U1,~,~] = svd(T1,"econ");
[U2,~,~] = svd(T2,"econ");
[U3,~,~] = svd(T3,"econ");


%% PART 2: Compute the core tensor

C = mode_n_product(T, U1', 1);
C = mode_n_product(C, U2', 2);
C = mode_n_product(C, U3', 3);

end