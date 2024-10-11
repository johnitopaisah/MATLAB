% Audio Signal Processing with Echo Cancellation
%
% Author: [Your Name]
% Date: [Date]
%
% Description:
% This MATLAB script performs echo cancellation on an audio signal stored
% in the file 'Pa11.wav'. It reads the audio signal, plots the original
% signal in the time and frequency domains, computes the Power Spectrum
% Density (PSD) of the signal, applies echo cancellation using a filter
% designed to remove echoes from the signal, and finally plots and plays
% the filtered signal.
%
% Usage:
% 1. Ensure the audio file 'Pa11.wav' is placed in the same directory as
%    this script.
% 2. Run the script. It will clear the command window, remove all variables
%    from the workspace, and close all open figures.
% 3. The original audio signal will be played first, followed by a pause
%    to allow for a moment of silence. Then, the processed (filtered)
%    audio signal will be played.
% 4. Listen to both signals to observe the difference in audio quality
%    before and after echo cancellation.
%
% Implementation Details:
% - The script reads the audio signal from the file 'Pa11.wav' and stores
%   it in the variable 'audio_signal', along with the sampling frequency
%   in 'sampling_frequency'.
% - It then computes the FFT of the audio signal to analyze its frequency
%   content and plots the original signal in both time and frequency
%   domains.
% - The script calculates the Power Spectrum Density (PSD) of the signal
%   and performs echo cancellation by designing a filter to remove echoes.
% - After applying the filter to the original signal, the filtered signal
%   is plotted in the time domain and played to hear the effect of echo
%   cancellation.
%
% Files:
% - `Pa11.wav`: Audio file used for demonstration.
%
% Dependencies:
% - MATLAB Signal Processing Toolbox (for FFT and inverse FFT).
%
% Note:
% - Ensure that the audio file 'Pa11.wav' is present in the same directory
%   as the script.
%


% Clear the command window, remove all variables from the workspace, and
% close all open figures
clc;
close all;

% Read the audio signal from the file 'Pa11.wav' and store the signal data
% in 'audio_signal' and the sampling frequency in 'sampling_frequency'.
[audio_signal, sampling_frequency] = audioread('Pa11.wav');

% Establish the signal length or number of samples
number_of_samples = length(audio_signal);

% Establish the signal sampling period
sampling_period = 1/sampling_frequency;

% Time vector for the audio signal
time_vector = 0:sampling_period:(number_of_samples-1)*sampling_period;

% Frequency range of the audio signal using the Nyquist frequency and the
% number of samples
nyquist_frequency = sampling_frequency/2;
frequency_interval = sampling_frequency/number_of_samples;
frequency_vector = 0:frequency_interval:(number_of_samples-1)*frequency_interval;

% Play the original signal
sound(audio_signal, sampling_frequency);
pause(length(audio_signal)/sampling_frequency); % Pause to allow for a moment of silence

% Plot the original audio signal in the time domain waveform
figure(1)
plot(time_vector, audio_signal)
title('Time Domain Waveform of Original Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Compute the FFT of the audio signal
fft_result = fft(audio_signal);

% Compute the magnitude of the FFT and shift the frequencies to center
% around zero
magnitude_spectrum = fftshift(abs(fft_result));

% Plot the magnitude spectrum of the original audio signal in the frequency
% domain
figure(2)
plot(frequency_vector, magnitude_spectrum)
title('Spectrum of the Original Audio Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Compute the Power Spectrum Density (PSD) of the signal
power_spectrum_density = magnitude_spectrum.^2;

% Plot the PSD
figure(3)
plot(frequency_vector, power_spectrum_density)
title('Power Spectrum Density (PSD) of the Signal');
xlabel('Frequency (Hz)')
ylabel('Power')

% Compute the inverse Fourier transform of the PSD to get the signal in
% time domain
time_domain_signal = ifft(power_spectrum_density);

% Establish a filter to remove echoes from the signal
[max_value_echo, index_echo] = max(time_domain_signal(1000:5000));
time_corresponding_echo = time_vector(index_echo-1);
start_index_filter = 1000 + index_echo - 1;
[max_value_original, index_original] = max(time_domain_signal(1:1000));
scale_factor = max_value_original/max_value_echo;
delta = (scale_factor*scale_factor) - 4;
filter_coefficient_rho = (scale_factor - sqrt(delta))/2;
filter_coefficient_array = zeros(1, start_index_filter + 1);
filter_coefficient_array(1) = 1;
filter_coefficient_array(end) = filter_coefficient_rho;

% Apply the filter to the original signal to get the filtered signal
filtered_signal = filter(1, filter_coefficient_array, audio_signal);

% Plot the filtered signal in the time domain
figure(4)
plot(time_vector, filtered_signal)
title('Filtered Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

% Play the filtered signal
sound(filtered_signal, sampling_frequency);
