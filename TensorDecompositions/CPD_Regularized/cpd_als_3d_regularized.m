function [B1, B2, B3, B4, c, output] = cpd_als_3d_regularized(T, M, R, options, varargin)
% CPD_ALS_3D Computes the canonical polyadic decomposition (CPD) of a
% 3rd-order tensor using the alternative least squares algorithm (ALS) 
% (assignment solution).
%
%INPUT:
%   T (I_1 x I_2 x I_3): 3-D tensor data tensor
%   R (1 x 1): number of CPD components
%   options (struct, optional) : optimization options containing:
%          - th_relerr (1 x 1): relative error threshold
%          - maxiter   (1 x 1): max number of iterations
%   variable extra inputs (1 x 3 cell array OR 3 arrays, OPTIONAL): 
%           initialization for the factor matrices
%
%OUTPUT:
%   B1 (I_1 x R): mode-1 factor matrix with normalized columns 
%   B2 (I_2 x R): mode-2 factor matrix with normalized columns 
%   B3 (I_3 x R): mode-3 factor matrix with normalized columns
%   c (1 x R):  vector containing the respective weights for the CPD
%               components
%   output (struct) : optimization options containing:
%          - numiter (1 x 1): the number of iterations the algorithm ran for
%          - relerr (1 x numiter): relative error achieved, defined as Frobenius norm of 
%           the residual of the decomposition OVER Frobenius norm of the 
%           original tensor.
%
% Remarks: 
%   The order in which the factor matrices is updated is fixed: mode-1, 
%   then mode-2 and, lastly, mode-3. 
%   
% Authors: Borbala Hunyadi (b.hunyadi@tudelft.nl)
%          Sofia-Eirini Kotti (s.e.kotti@tudelft.nl)

% Get algorithm parameters
% Number of iterations
maxiter = options.maxiter;
% Threshold for relative error
th_relerr = options.th_relerr;
% regularization parameter rho
rho = options.rho;
%
% Check if the initialization for the factor matrices was given
init = [];
if ~isempty(varargin)
    if length(varargin) == 1    % Given as cell
        init = varargin{:}; 
    else                        % Given as matrices 
        init = varargin;
    end
end

% Initialize the three factor matrices 
if isempty(init)    % Randomly if not initialization was given
    B1 = randn(size(T, 1), R);
    B2 = randn(size(T, 2), R);
    B3 = randn(size(T, 3), R);
    B4 = randn(R, R);
else                % Otherwise use the given initialization
    assert(all(cellfun(@(x) size(x,2), init) == R), ...
        "The given initialization has a different number of rank-1 components than the given R.")
    B1 = init{1};
    B2 = init{2};
    B3 = init{3};
    B4 = init{4};
end


relerr = Inf;
c = sqrt(sum(B1.*B1,1)).*sqrt(sum(B2.*B2,1)).*sqrt(sum(B3.*B3,1));

B1 = 1./sqrt(sum(B1.*B1,1)).*B1;
B2 = 1./sqrt(sum(B2.*B2,1)).*B2;
B3 = 1./sqrt(sum(B3.*B3,1)).*B3;
B4 = 1./sqrt(sum(B4.*B4,1)).*B4;


% Obtain the three tensor unfoldings
T1 = mode_n_matricization(T, 1);
T2 = mode_n_matricization(T, 2);
T3 = mode_n_matricization(T, 3);


% ALS iterations
for idxiter = 1:maxiter

    % Mode 1 (do not forget to normalize the columns!)
    B1 = T1 * khatri_rao(B3, B2) * pinv((B3'*B3).*(B2'*B2));
    B1 = 1./sqrt(sum(B1.*B1,1)).*B1;

    % Mode 2 (do not forget to normalize the columns!)
    B2 = (T2 * khatri_rao(B3, B1) + rho*(M*B2*B4' + M'*B2*B4)) * pinv((B3'*B3).*(B1'*B1) + rho*(B4'*(B2')*B2*B4 + B4*(B2')*B2));
    % B2 = (T2 * khatri_rao(B3, B1) + rho*var1) * pinv((B3'*B3).*(B1'*B1) + rho*var2);
    B2 = 1./sqrt(sum(B2.*B2,1)).*B2;

    % Mode 3 (do not forget to normalize the columns!)
    B3 = T3 * khatri_rao(B2, B1) * pinv((B2'*B2).*(B1'*B1));
    c = sqrt(sum(B3.*B3,1));
    B3 = 1./sqrt(sum(B3.*B3,1)).*B3;

    % B4 = (B2' * M * B2) * pinv((B2' * B2) * B4 * (B2' * B2));
    B4 = inv(B2' * B2)*B2' * M * B2 * inv(B2' * B2);
    B4 = 1./sqrt(sum(B4.*B4,1)).*B4;

    % Compute the current estimate of the mode-3 unfolding of the tensor
    T3_hat = B3 * diag(c) * khatri_rao(B2, B1)';

    % Calculate the relative error between the estimate and the true
    % mode-3 unfolding of the tensor
    % Update the following line
    relerr(idxiter) = (1/2) * (norm(T - reshape(B1 * khatri_rao(B3, B2)', size(T)), "fro")^2 + (rho/2) * norm(M - B2 * B4 * B2', "fro")^2);

    % Check stopping criterion on relative error
    if relerr(idxiter) < th_relerr 
        break;
    end

end

% Add the relative error and the number of iterations to output structure
output.relerr = relerr;
output.numiter = idxiter;
% 
% Warning if maximum number of iterations was reached
if idxiter == options.maxiter
    warning(['The ALS algorithm reached the maximum number of ' num2str(options.maxiter) ' iterations. Error: ' num2str(relerr(end))])
end
