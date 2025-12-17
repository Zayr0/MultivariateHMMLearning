function [pcaData, coeff, latent, meanData] = run_pca(data, varargin)
% RUN_PCA Perform PCA on channels × samples data.
%
% INPUT:
%   data : channels × samples
%
% OPTIONAL:
%   'nPC'       : number of principal components to keep (default: all)
%   'varThresh' : percent variance to keep (0–1), overrides 
%   'whiten'    : true/false (default: false)nPC
%
% OUTPUT:
%   pcaData : projected data (nPC × samples)
%   coeff   : PCA loading matrix (channels × nPC)
%   latent  : eigenvalues
%   meanData: channel mean removed before PCA
%
% Example:
%   [pc, W] = run_pca(eeg, 'nPC', 20, 'whiten', true);

    p = inputParser;
    addParameter(p, 'nPC', size(data,1));
    addParameter(p, 'varThresh', []);
    addParameter(p, 'whiten', false);
    parse(p, varargin{:});

    nPC        = p.Results.nPC;
    varThresh  = p.Results.varThresh;
    whitenFlag = p.Results.whiten;

    % Center the data
    meanData = mean(data, 2);
    X = data - meanData;

    % Covariance PCA
    C = cov(X.');
    [V, D] = eig(C);
    [latent, idx] = sort(diag(D), 'descend');
    coeff = V(:, idx);

    % Determine number of PCs
    if ~isempty(varThresh)
        cumVar = cumsum(latent) / sum(latent);
        nPC = find(cumVar >= varThresh, 1, 'first');
    end

    coeff = coeff(:, 1:nPC);
    latent = latent(1:nPC);

    % Project data
    pcaData = coeff.' * X;

    % Optional whitening
    if whitenFlag
        pcaData = diag(1 ./ sqrt(latent)) * pcaData;
    end
end
