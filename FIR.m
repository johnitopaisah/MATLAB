%--------------------------------------------------------------------------
%    Audio Signal Filtering and Analysis Script
%
%    Part 2: FIR filter as moving average filter
%
%    STUDENT_1 NAME: John Itopa ISAH - 63362
%    STUDENT_2 NAME: Blesso Danny Jeba Stephenen - 63291
%
%
% This script performs filtering and spectral analysis on an audio signal
% using a moving average filter (Gaussian window) with different window
% sizes (N = 3, 11, 31). It also visualizes the time and frequency domain
% analysis of the original and filtered signals.
%
% Usage:
%   1. Ensure the audio file 'Don_Giovanni_2.wav' is located in the same
%      directory as this script.
%   2. Run the script in MATLAB.
%
% Implementation Details:
%
% 1. Load the audio signal from the file 'Don_Giovanni_2.wav' and extract
%    the sampling frequency and signal length.
%
% 2. Define parameters for time and frequency vectors based on the sampling
%    frequency and signal length.
%
% 3. Iterate over different window sizes (N = 3, 11, 31):
%    - Compute the Gaussian window function for each window size.
%    - Apply the digital filter to the input audio signal using the computed
%      window function.
%    - Play the original and filtered audio signals sequentially with a
%      one-second delay between them.
%
% 4. Compute and plot the time domain analysis:
%    - Plot the original and filtered audio signals in the time domain.
%
% 5. Compute the FFT (Fast Fourier Transform) of the original and filtered
%    audio signals to obtain their magnitude spectra.
%
% 6. Plot the frequency domain analysis:
%    - Plot the magnitude spectrum of the original and filtered signals in
%      the frequency domain.
%
% 7. Display the plots for each window size (N = 3, 11, 31) in separate
%    figures.
%
% Dependencies:
%   - MATLAB Signal Processing Toolbox
%
% Note:
%   - Ensure the audio file 'Don_Giovanni_2.wav' is present in the same
%     directory as this script for proper execution.
%
%--------------------------------------------------------------------------

clc; % Clear the Command Window
close all; % Close all open figures

% Load audio signal and extract parameters
[audioSignal, samplingFreq] = audioread('Don_Giovanni_2.wav');
signalLength = length(audioSignal);
samplingPeriod = 1 / samplingFreq;
timeVector = 0:samplingPeriod:(signalLength-1)*samplingPeriod;
frequencyInterval = samplingFreq / signalLength;
frequencyVector = 0:frequencyInterval:(signalLength-1)*frequencyInterval;

% Play original audio signal with delay
sound(audioSignal, samplingFreq);
pause(signalLength / samplingFreq + 1); % Delay for the length of the audio + 1 second
    

% Define different window sizes (N values)
N_values = [3, 11, 31];

% Iterate over each window size
for i = 1:length(N_values)
    % Current window size
    windowSize = N_values(i);
    
    % Compute Gaussian window function
    windowFunction = gausswin(2*windowSize+1) ./ windowSize;

    % Apply the digital filter to the audio signal
    filteredSignal = filter(windowFunction, 1, audioSignal);
    
    % Play filtered audio signal with delay
    sound(filteredSignal, samplingFreq);
    pause(signalLength / samplingFreq + 1); % Delay for the length of the audio + 1 second
    
    % Compute FFT of the input audio signal and filtered signal
    fftMagnitude = abs(fft(audioSignal));
    filteredFFTMagnitude = abs(fft(filteredSignal));
    
    % Plot time domain analysis
    figure;
    plot(timeVector, audioSignal, 'r');
    hold on;
    plot(timeVector, filteredSignal, 'b');
    title(sprintf('Time Domain Analysis N=%d', windowSize));
    xlabel('Time in seconds');
    ylabel('Amplitude');
    legend('Before Filtering Time', 'After Filtering Time');
    grid on;
    
    % Plot frequency domain analysis
    figure;
    semilogy(frequencyVector, fftMagnitude, 'y');
    hold on;
    semilogy(frequencyVector, filteredFFTMagnitude, 'g');
    title(sprintf('Spectrum of Signal N=%d', windowSize));
    xlabel('Frequency in Hz');
    ylabel('Amplitude in log scale');
    legend('Before Filtering Freq', 'After Filtering Freq');
    grid on;
end
