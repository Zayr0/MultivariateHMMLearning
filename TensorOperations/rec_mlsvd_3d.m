function T_rec = rec_mlsvd_3d(Ct,U1t,U2t,U3t)

% REC_MLSVD_3D reconstructs a 3-D tensor from its MLSVD representation.
%
% INPUT:
%   Ct (3D array)    : MLSVD core tensor (truncated)
%   U1t (matrix)     : MLSVD factor matrix of the first mode (truncated)
%   U2t (matrix)     : MLSVD factor matrix of the second mode (truncated)
%   U3t (matrix)     : MLSVD factor matrix of the third mode (truncated)
% 
% OUTPUT:
%   T_rec (3D array) : Tensor reconstructed from the truncated MLSVD
% expression. It gives an approximation of the original tensor (T).


T_rec = mode_n_product(Ct, U3t, 3);
T_rec = mode_n_product(T_rec, U2t, 2);
T_rec = mode_n_product(T_rec, U1t, 1);
end