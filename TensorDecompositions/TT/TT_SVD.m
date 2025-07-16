function [tt_cores,rel_error] = TT_SVD(tensor, epsilon)
%TT_SVD(tensor,epsilon) 
%   Algorithm that decomposes a given N-th order tensor into tensor train format
%INPUT:
%   tensor (N-dimensional double):  N-th order tensor
%   epsilon (double):               error-bound for approximation error of TT decomposition
%OUTPUT:
%   tt_cores (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%   error (double):                 actual approximation error of TT
%                                   decomposition
n = ndims(tensor);
% every cell(1,n) contains one tt-core, which is a r_i x n_i x r_i+1 double
% do not change the datastructure! Otherwise, the testfunction is not
% guaranteed to work.
% The cell(1,n+1) contains the location of the core at the moment. 
% If the norm is stored in the 10th core, the value in the cell(1,d+1) should be 10. 
% If unsure, have a look at the variable tt_dog
tt_cores = cell(1,n+1);  
%------ Implement your code below ------

rel_error = 0;

Fnorm2 = reshape(tensor,1,[])*reshape(tensor,[],1);

Lo = tensor;
R = 1;

for i = 1:(n-1)
    Lo = reshape(Lo, R*size(tensor,i), []);
    [U, S, V] = svd(Lo, 'econ');
    
    error = 0;
    tID = size(S, 1);
    for j = size(S,1):-1:2
        temp_error = diag(S(j:end,j:end))' * diag(S(j:end,j:end));
        if temp_error <= (epsilon^2) * Fnorm2 / (n-1)
            tID = j-1;
            error = temp_error;
            break
        end
    end
    S = S(1:tID,1:tID);
    U = U(:,1:tID);
    V = V(:,1:tID);

    rel_error = rel_error + error/Fnorm2;
    Lo = S * V';
    R = size(S, 1);

    if i == 1
        tt_cores{i} = U; 
    else
        tt_cores{i} = reshape(U,[],size(tensor, i),R);
    end
end
tt_cores{n} = S * V';
tt_cores{n+1} = n;
rel_error = sqrt(rel_error);
end