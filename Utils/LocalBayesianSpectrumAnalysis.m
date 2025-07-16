function [tVec, ampVec] = LocalBayesianSpectrumAnalysis(x, fs, winLength, overlap)
% localBayesianSpectrum - Estimates local amplitude via Bayesian-inspired spectral analysis
% 
% Syntax:
%   [tVec, ampVec] = localBayesianSpectrum(x, fs, winLength, overlap)
%
% Inputs:
%   x         - Input time series (1D vector)
%   fs        - Sampling frequency (Hz)
%   winLength - Window length (samples)
%   overlap   - Overlap between windows (samples)
%
% Outputs:
%   tVec      - Time vector (center of each window)
%   ampVec    - Estimated amplitude trajectory

% Default to Hamming window
w = hamming(winLength);

% Segment signal
step = winLength - overlap;
nWins = floor((length(x) - winLength) / step) + 1;

ampVec = zeros(nWins, 1);
tVec = zeros(nWins, 1);

for i = 1:nWins
    idxStart = (i-1)*step + 1;
    idxEnd = idxStart + winLength - 1;
    segment = x(idxStart:idxEnd) .* w;

    % FFT
    X = fft(segment);
    X = X(1:floor(winLength/2));
    f = (0:length(X)-1)*(fs/winLength);

    % Compute power spectrum
    P = abs(X).^2;

    % Bayesian-inspired: Normalize to probability distribution
    P_norm = P / sum(P + eps);

    % Weighted amplitude: Emphasize peak frequency content
    amp = sum(P .* P_norm);

    % Store results
    ampVec(i) = sqrt(amp);  % sqrt to get back to amplitude units
    tVec(i) = (idxStart + winLength/2 - 1) / fs;
end

end
