function [freqs, timeVec] = localFrequencyAnalysis(s, fs, windowSize, overlap, method, smoothing)
% localFrequencyAnalysis - Computes dominant local frequency over time
% Syntax:
%   [freqs, timeVec] = localFrequency(s, fs, windowSize, overlap, method)
%
% Inputs:
%   s          - Input signal (1D array)
%   fs         - Sampling frequency (Hz)
%   windowSize - Window size in seconds (e.g., 0.1)
%   overlap    - Overlap between windows in seconds (e.g., 0.05)
%   method     - 'average' or 'dominant'
%
% Outputs:
%   freqs   - Local frequency in each window (Hz)
%   timeVec - Time vector corresponding to each window center

    if size(s, 1) > 1
        s = s(:); % Ensure column vector
    end

    if nargin < 6
        smoothing = 0; % No smoothing by default
    end

    % Convert to sample units
    winSamples = round(windowSize * fs);
    overlapSamples = round(overlap * fs);

    % Compute STFT
    [S, F, T] = stft(s, fs, ...
        'Window', hann(winSamples, 'periodic'), ...
        'OverlapLength', overlapSamples, ...
        'FFTLength', 2^nextpow2(winSamples));

    % Power spectrum
    P = abs(S).^2;

    switch lower(method)
        case 'dominant'
            % Get frequency of max power per window
            [~, maxIdx] = max(P, [], 1);
            freqs = abs(F(maxIdx));

        case 'average'
            % Normalize power to use as weights
            P_sum = sum(P, 1);
            P_sum(P_sum == 0) = eps;  % Avoid divide-by-zero
            P_norm = P ./ P_sum;

            % Compute weighted average frequency
            freqs = sum(P_norm .* F, 1);

        otherwise
            error('Invalid method. Use ''dominant'' or ''average''.');
    end

    if smoothing > 1
        freqs = smoothdata(freqs, 'movmean', smoothing);
    end

    timeVec = T;
end
