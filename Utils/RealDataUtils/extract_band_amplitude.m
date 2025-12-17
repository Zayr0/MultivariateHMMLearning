function [ampEnv, envTime] = extract_band_amplitude(data, fs, band, varargin)
% EXTRACT_BAND_AMPLITUDE 
%   Computes amplitude envelope (Hilbert-based) for a specific frequency band.
%
% INPUTS:
%   data : channels × samples matrix (EEG/EOG data)
%   fs   : sampling rate (Hz)
%   band : [lowFreq highFreq] frequency band
%
% OPTIONAL PARAMETERS:
%   'ds'        : downsample factor (default = 1 = no downsampling)
%   'filtOrder' : FIR filter order (default = fs*2)
%
% OUTPUTS:
%   ampEnv  : amplitude envelope (channels × samples/downsample)
%   envTime : time vector for amplitude envelope
%
% Example:
%   [env, t] = extract_band_amplitude(eeg, 500, [8 12], 'ds', 4);
%
% ---------------------------------------------------------------

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'ds', 1, @(x) isnumeric(x) && x >= 1);
    addParameter(p, 'filtOrder', fs*2, @(x) isnumeric(x) && x > 0);
    parse(p, varargin{:});
    ds = p.Results.ds;
    filtOrder = p.Results.filtOrder;

    % Ensure row-major for channels
    if size(data,1) > size(data,2)
        warning('Data appears to be samples × channels; transposing.');
        data = data.';
    end

    % --- 1. Band-pass filter ------------------------------------------------
    filt = designfilt('bandpassfir', ...
                      'FilterOrder', filtOrder, ...
                      'CutoffFrequency1', band(1), ...
                      'CutoffFrequency2', band(2), ...
                      'SampleRate', fs);

    % Zero-phase filtering
    filtered = filtfilt(filt, data.').';  % preserve channels × samples

    % --- 2. Hilbert transform to compute analytic signal ---------------------
    analyticSig = hilbert(filtered.').';
    ampEnvFull  = abs(analyticSig);

    % --- 3. Downsample envelope if requested --------------------------------
    if ds > 1
        ampEnv = ampEnvFull(:, 1:ds:end);
    else
        ampEnv = ampEnvFull;
    end

    % Time vector
    envTime = (0:size(ampEnv,2)-1) * (1/fs) * ds;
end
