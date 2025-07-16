function activitySignal = localSpectralActivity(signal, fs, windowSize, overlap, freqBand)
% localSpectralActivity - Computes a signal scaled by local frequency activity
%
% Inputs:
%   signal     - Input time series (1D vector)
%   fs         - Sampling frequency (Hz)
%   windowSize - Window size in samples for local analysis (e.g., 256)
%   overlap    - Overlap in samples (e.g., 128)
%   freqBand   - Frequency band of interest [fLow, fHigh] in Hz
%
% Output:
%   activitySignal - Amplitude modulated signal reflecting local spectral activity

    % Ensure column vector
    signal = signal(:);

    % Define STFT parameters
    nfft = max(256, 2^nextpow2(windowSize));
    win = hamming(windowSize);

    % Compute spectrogram (magnitude squared)
    [S, F, T] = spectrogram(signal, win, overlap, nfft, fs, 'yaxis');
    P = abs(S).^2;  % Power spectrum

    % Frequency band indices
    bandIdx = F >= freqBand(1) & F <= freqBand(2);

    % Band power over time
    bandPower = sqrt(sum(P(bandIdx, :), 1));  % root-bandpower for amplitude scaling

    % Interpolate bandPower to match original signal length
    timeIdx = linspace(1, length(signal), length(bandPower));
    bandPowerInterp = interp1(timeIdx, bandPower, 1:length(signal), 'linear', 'extrap');

    % Scale original signal by local band activity
    activitySignal = signal .* bandPowerInterp(:);
end
