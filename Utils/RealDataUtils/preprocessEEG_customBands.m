function [phaseData, ampEnvData, dAmpEnvData] = preprocessEEG_customBands(eegData, fs, freqBands, varargin)
% preprocessEEG_customBands: Preprocess EEG data with flexible options
%
% Inputs:
%   eegData   : matrix of EEG signals [channels x timepoints]
%   fs        : sampling rate in Hz
%   freqBands : Nx2 matrix of frequency ranges (Hz), e.g., [1 3; 7 9; 13 15]
%
% Optional Name-Value Inputs:
%   'FilterOrder' : FIR filter order (default = 3000)
%   'DT'          : time step in seconds for derivative (default = 0.002)
%
% Outputs:
%   phaseData    : cell array containing phase of filtered signals for each band
%   ampEnvData   : cell array containing amplitude envelope for each band
%   dAmpEnvData  : cell array containing time derivative of amplitude envelope

%% Parse optional inputs
p = inputParser;
addParameter(p, 'FilterOrder', 3000, @(x) isnumeric(x) && x > 0);
addParameter(p, 'DT', 0.002, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});

filterOrder = p.Results.FilterOrder;
dt = p.Results.DT;

%% Convert freqBands if it is a vector
if isvector(freqBands)
    % Make consecutive pairs
    freqBands = [freqBands(1:end-1)' freqBands(2:end)'];
end

%% Initialize
nBands = size(freqBands,1);
nChannels = size(eegData,1);

phaseData = cell(nBands,1);
ampEnvData = cell(nBands,1);
dAmpEnvData = cell(nBands,1);

diffFactor = dt * fs;

%% Loop over frequency bands
for b = 1:nBands
    fLow = freqBands(b,1);
    fHigh = freqBands(b,2);

    % Design FIR bandpass filter
    bFIR = fir1(filterOrder, [fLow fHigh]/(fs/2), 'bandpass');

    % Zero-phase filtering
    filteredData = filtfilt(bFIR, 1, double(eegData') )';

    % Hilbert transform
    analyticSignal = hilbert(filteredData')';

    % Phase and amplitude envelope
    phaseData{b} = angle(analyticSignal);
    ampEnvData{b} = abs(analyticSignal);

    % Time derivative of amplitude envelope
    dAmpEnvData{b} = [diff(ampEnvData{b},1,2) zeros(nChannels,1)] / diffFactor;
end

end
