function [icaData, A, W] = run_ica(data, varargin)
% RUN_ICA Perform ICA on channels × samples EEG/EOG data.
%
% INPUT:
%   data : channels × samples
%
% OPTIONAL:
%   'nIC' : number of independent components to extract (default: all)
%
% OUTPUT:
%   icaData : independent components (nIC × samples)
%   A       : mixing matrix   (channels × nIC)
%   W       : unmixing matrix (nIC × channels)
%
% Example:
%   [ic, A, W] = run_ica(eeg, 'nIC', 30);

    p = inputParser;
    addParameter(p, 'nIC', size(data,1));
    parse(p, varargin{:});
    nIC = p.Results.nIC;

    % Center
    meanData = mean(data, 2);
    X = data - meanData;

    % Run FastICA with symmetric approach
    [icaData, A, W] = fastica(X, ...
                              'numOfIC', nIC, ...
                              'approach', 'symm', ...
                              'verbose', 'off');
end
