%  Dual-Tone Multi-Frequency (DTMF) Signal Detection
%
%   STUDENT_1 -> John Itopa ISAH - 63362
%   STUDENT_2 -> Blesso Danny Jeba Stephen - 63291
%
% Description:
% This MATLAB script detects Dual-Tone Multi-Frequency (DTMF) signals in an
% audio file. It reads the input audio signal from a WAV file, performs
% spectral analysis to identify DTMF tones, applies Chebyshev Type I filters
% to isolate relevant frequency bands, and maps the detected tones to the
% corresponding dialed numbers (0-9, *, #, A-D). The script provides
% functionality to plot the original signal in the time and frequency
% domains, as well as to play the original and filtered signals.
%
% Usage:
% 1. Ensure the desired input audio files are saved in WAV format and placed
%    in the same directory as this script.
% 2. Modify the filename in the script to match the name of your input audio
%    file, if different from the default ('0123456789.wav').
% 3. Run the script in MATLAB.
% 4. Listen to the output audio to hear the detected numbers corresponding
%    to the dialed DTMF tones.
%
% Implementation Details:
% - The script reads the input audio signal from a WAV file and calculates
%   relevant sampling parameters such as sampling frequency and signal length.
% - It transforms the audio signal into the frequency domain using the Fast
%   Fourier Transform (FFT) and performs spectral analysis on each window
%   of the signal.
% - Chebyshev Type I filters are designed and applied to isolate the relevant
%   frequency bands corresponding to DTMF tones.
% - Detected frequency pairs are mapped to the corresponding dialed numbers
%   based on predefined frequency ranges with a 1.5% frequency tolerance.
%
% Files:
% - `DTMF_SignalDetection.m`: MATLAB script for DTMF signal detection.
% - Sample audio files (e.g., 'numero3.wav') for demonstration.
%
% Dependencies:
% - MATLAB Signal Processing Toolbox (for FFT, filtering, and spectral analysis).
%
% Note:
% - Ensure that the input audio file is present in the same directory as
%   the script and is in WAV format.



% Read the audio file
[audio_signal, sampling_frequency] = audioread('0123456789.wav');
%[audio_signal, sampling_frequency] = audioread('numero2.wav');
%[audio_signal, sampling_frequency] = audioread('numero3.wav');
%[audio_signal, sampling_frequency] = audioread('numero4.wav');
%sound(audio_signal)
% Compute the sampling period Ts
sampling_period = 1 / sampling_frequency;
% Compute the signal length or number of samples in the signal
number_of_samples = length(audio_signal);
% Create time vector for the audio signal.
time_vector = 0:sampling_period:(number_of_samples - 1) * sampling_period;

% Plot the original audio signal of in time domain waveform
figure(1)
plot(time_vector, audio_signal)
xlabel('Time in seconds')
ylabel('Sampled Signal')
title('Given Signal in Time Domain')

% Compute the frequency range of the audio signal using the Nyquist frequency and the number of samples
% f = sampling_frequency / 2; % Establish the Nyquist frequency 
% f_interval = sampling_frequency / number_of_samples; % Interval between frequency samples
% f = 0:f_interval:(number_of_samples - 1) * f_interval; % Create frequency vector using the interval and number of samples

%% Transform the signal into the frequency domain with FFT and store the magnitude of the FFT in X_mag.
% Y = fft(audio_signal); % Compute the FFT of the audio signal
% Y_mag = fftshift(abs(Y)); % Compute the magnitude of the FFT and shift the frequencies to center around zero

f = sampling_frequency / 2; % Establish the Nyquist frequency
% 256 is chosen for an optimum frequency resolution for the spectrum 
M = 256;
% Create a linearly spaced frequency vector from 0 Hz to sampling frequency
frequency_vector = linspace(0, sampling_frequency, M + 1);
% Remove the last element from the frequency vector and discard the Nyquist frequency point from the vector.
frequency_vector = frequency_vector(1:end - 1);
% Convert the frequencies from the range of 0 Hz to the Nyquist frequency to the normalized range of 0 to 1, where 1 represents the Nyquist frequency.
normalized_frequency = frequency_vector / sampling_frequency;
% Shift the frequencies in the vector frequency_vector to be centered around zero
frequency_vector = frequency_vector - f;
% Compute the FFT of the signal with resolution M to examine the magnitude and phase information in the resulting spectrum
fft_signal = fft(audio_signal, M);

figure(2) 
plot(normalized_frequency, abs(fft_signal)) % Plot the magnitude spectrum of the original audio signal in frequency domain
title('Spectrum of the original signal before filtering');
xlabel('Frequency (Hz)');
ylabel('Spectrum of the Signal')
grid on;

%% Spectrum analysis of each window in 65ms
window_length = 65e-3;
% Determine the number of samples for each window
samples_per_window = round(window_length / sampling_period);

% Reshaping the signal in intervals of the window size
% N = length(audio_signal);
f_interval = 1 / number_of_samples; % Interval between frequency samples
frequency_vector = 0:f_interval:(number_of_samples - 1) * f_interval;% Create frequency vector containing the frequencies corresponding to each frequency 
% sample point in the Fourier transform

% Reshaping the signal in intervals of the window size
number_of_zeros = ceil(number_of_samples / samples_per_window) * samples_per_window - number_of_samples; % Determine the number of zeros that need to be added to the signal
padded_signal = [audio_signal.' zeros(1, number_of_zeros)];% Concatenate the original signal audio_signal with a vector of zeros of length number_of_zeros
reshaped_signal = reshape(padded_signal', samples_per_window, [])'; % Reshape the signal by transposing padded_signal
filtered_signal = zeros(size(reshaped_signal)); % Initialize a matrix of zeros with the same size as reshaped_signal

time_vector_padded = [time_vector  zeros(1, number_of_zeros)];% Concatenate the time vector time_vector with a vector of zeros of length number_of_zeros
time_vector_reshaped = reshape(time_vector_padded', samples_per_window, [])'; % Restore the original orientation and has the same dimensions as the reshaped signal matrix reshaped_signal

%% Minimum power initialization
minimum_power_threshold = 20;

% Initialization of the vector
selected_windows = zeros(size(reshaped_signal,1), 1); % Create a placeholder for storing results related to the processing of the signal segments

% Iterate over the rows of the reshaped_signal matrix
for k = 1:size(reshaped_signal, 1)
    
    power = sum(reshaped_signal(k, :).^2) / samples_per_window; % Power of the signal
    power_dB = (10 * log10(power / 10.^-3)); % Compute power of the signal in dB
    
    selected_windows(k) = power_dB > minimum_power_threshold; % Select all windows with power greater 20dB
    
end
% Compute the number of windows or signal segments that have power levels greater than the threshold of 20dB
number_of_selected_windows = sum(selected_windows);

power = zeros(size(reshaped_signal, 1), 1); % Placeholder used to store the calculated power of each signal segment
power_dB = zeros(size(reshaped_signal, 1), 1); % Used to store the power values converted to decibels (dB)
frequency_pairs = zeros(number_of_selected_windows, 2); % Used to store the frequencies obtained from the signal segments that have power levels above the specified threshold

%% Parameters for cheby1() filter function
filter_order = 7;
stopband_attenuation = 10;
frequency_range_factor_1 = 0.55;
frequency_range_factor_2 = 0.93;

% Iterating over rows of reshaped_signal
for k = 1:size(reshaped_signal, 1)
    
    if  k == 1      %
        % Calculate the Fast Fourier Transform (FFT) of the k-th row of reshaped_signal using N points.
        FFT_signal_1 = fft(reshaped_signal(k, :), number_of_samples);
        % Compute the magnitude of the FFT result
        FFT_signal_1_magnitude = abs(FFT_signal_1);
        % Find the peaks in the magnitude spectrum of the positive frequencies
        [peaks_1, locations_1] = findpeaks(FFT_signal_1_magnitude(1:(number_of_samples + 1) / 2), 'SortStr', 'descend');
        % Assign the frequency corresponding to the highest peak to the first column of frequency_pairs.
        frequency_pairs(k, 1) = frequency_vector(locations_1(1)); % Get the location of the first frequency
        normalized_frequency_1 = frequency_pairs(k, 1); % Stores the frequency value

        % Design a Chebyshev Type I filter with the specified parameters based on the frequency range [frequency_range_factor_1 * normalized_frequency_1 frequency_range_factor_2 * normalized_frequency_1]
        [b, a] = cheby1(filter_order, stopband_attenuation, [frequency_range_factor_1 * normalized_frequency_1 frequency_range_factor_2 * normalized_frequency_1], 'stop'); 
        % Apply the designed filter b and a to the k-th row of reshaped_signal and store the filtered result
        filtered_signal(k, :) = filter(b, a, reshaped_signal(k, :));
        % Compute the FFT of the filtered signal
        FFT_signal_2 = fft(filtered_signal(k, :), number_of_samples);
        % Calculate the magnitude spectrum of the filtered signal.
        FFT_signal_2_magnitude = abs(FFT_signal_2);
        % Find the peaks in the magnitude spectrum of the positive frequencies 
        [peaks_2, locations_2] = findpeaks(FFT_signal_2_magnitude(1:(number_of_samples + 1) / 2), 'SortStr', 'descend');
        % Assign the frequency corresponding to the highest peak to the second column of frequency_pairs
        frequency_pairs(k, 2) = frequency_vector(locations_2(1));
        
    else
        % For iterations other than the first one
        power(k - 1) = sum(reshaped_signal(k - 1, :).^2) / samples_per_window; % Power of the signal
        power_dB(k - 1) = (10 * log10(power(k - 1) / 10.^-3)); % Power of the signal in DB
        
        power(k) = sum(reshaped_signal(k, :).^2) / samples_per_window; % Power of the signal
        power_dB(k) = (10 * log10(power(k) / 10.^-3)); % Power of the signal in DB
        
        if power_dB(k - 1) < minimum_power_threshold && power_dB(k) > minimum_power_threshold % To avoid taking duplicate frequencies
            
            FFT_signal_1 = fft(reshaped_signal(k, :), number_of_samples); % Calculate the FFT of the k-th row of reshaped_signal using N points.
            
            FFT_signal_1_magnitude = abs(FFT_signal_1); % Compute the magnitude of the FFT result
            % Find the peaks in the magnitude spectrum of the positive frequencies up to Nyquist frequency.
            [peaks_1, locations_1] = findpeaks(FFT_signal_1_magnitude(1:(number_of_samples + 1) / 2), 'SortStr', 'descend');
            % Assign the frequency corresponding to the highest peak to the first column of frequency_pairs.
            frequency_pairs(k, 1) = frequency_vector(locations_1(1));
   
            normalized_frequency_1 = frequency_pairs(k, 1); % Stores the frequency value in normalized_frequency_1.
            [b, a] = cheby1(filter_order, stopband_attenuation, [frequency_range_factor_1 * normalized_frequency_1 frequency_range_factor_2 * normalized_frequency_1], 'stop');% To design the filter
            
            filtered_signal(k, :) = filter(b, a, reshaped_signal(k, :)); % To apply the filter 
            FFT_signal_2 = fft(filtered_signal(k, :), number_of_samples); % Compute the FFT of the filtered signal
            FFT_signal_2_magnitude = abs(FFT_signal_2); % Calculate the magnitude spectrum of the filtered signal.
            % Find the peaks in the magnitude spectrum of the positive frequencies up to Nyquist frequency.
            [peaks_2, locations_2] = findpeaks(FFT_signal_2_magnitude(1:(number_of_samples + 1) / 2), 'SortStr', 'descend');
            frequency_pairs(k, 2) = frequency_vector(locations_2(1)); % Get the location of the second frequency
            
        end
    end
end

% The first and second column were sorted to make the first column to be
% the low frequency and second column the high frequency
sorted_frequency_pairs = reshape(nonzeros(frequency_pairs), [], 2) * sampling_frequency; % Obtain the frequencies in Hz
sorted_frequency_pairs = sort(sorted_frequency_pairs.', 'ascend'); % Transpose Tone3 matrix and sorts it in ascending order, resulting in a 2-column matrix for the frequency pair
sorted_frequency_pairs = sorted_frequency_pairs.'; % Restore the original arrangement of frequency pairs
detected_numbers = zeros(length(sorted_frequency_pairs), 1); % To store the calculated values (numbers dialed) in the subsequent code

%% To compare the frequencies in the signal to the DTMF, considering the 1.5% frequency tolerance

for k = 1:length(sorted_frequency_pairs)
    
    if ~isempty(sorted_frequency_pairs) && size(sorted_frequency_pairs, 2) >= 2 && sorted_frequency_pairs(k, 1) >= 686.545 && sorted_frequency_pairs(k, 1) <= 707.455 && sorted_frequency_pairs(k, 2) >= 1190.865 && sorted_frequency_pairs(k, 2) <= 1227.135
        detected_numbers(k) = 1;
        
    elseif  sorted_frequency_pairs(k, 1) >= 686.545 && sorted_frequency_pairs(k, 1) <= 707.455 && sorted_frequency_pairs(k, 2) >= 1315.96 && sorted_frequency_pairs(k, 2) <= 1356.04
        detected_numbers(k) = 2;
        
    elseif sorted_frequency_pairs(k, 1) >= 686.545 && sorted_frequency_pairs(k, 1) <= 707.455 && sorted_frequency_pairs(k, 2) >= 1454.845 && sorted_frequency_pairs(k, 2) <= 1499.155
        detected_numbers(k) = 3;
        
    elseif  sorted_frequency_pairs(k, 1) >= 686.545 && sorted_frequency_pairs(k, 1) <= 707.455 && sorted_frequency_pairs(k, 2) >= 1612.445 && sorted_frequency_pairs(k, 2) <= 1661.555
        detected_numbers(k) = 'A';
        
    elseif sorted_frequency_pairs(k, 1) >= 758.45 && sorted_frequency_pairs(k, 1) <= 781.55 && sorted_frequency_pairs(k, 2) >= 1190.865 && sorted_frequency_pairs(k, 2) <= 1227.135
        detected_numbers(k) = 4;
        
    elseif  sorted_frequency_pairs(k, 1) >= 758.45 && sorted_frequency_pairs(k, 1) <= 781.55  && sorted_frequency_pairs(k, 2) >= 1315.96 && sorted_frequency_pairs(k, 2) <= 1356.04
        detected_numbers(k) = 5;
        
    elseif sorted_frequency_pairs(k, 1) >= 758.45 && sorted_frequency_pairs(k, 1) <= 781.55  && sorted_frequency_pairs(k, 2) >= 1454.845 && sorted_frequency_pairs(k, 2) <= 1499.155
        detected_numbers(k)= 6;
        
    elseif  sorted_frequency_pairs(k, 1) >= 758.45 && sorted_frequency_pairs(k, 1) <= 781.55  && sorted_frequency_pairs(k, 2) >= 1612.445 && sorted_frequency_pairs(k, 2) <= 1661.555
        detected_numbers(k) = 'B';
        
    elseif sorted_frequency_pairs(k, 1) >= 839.22 && sorted_frequency_pairs(k, 1) <= 864.78 && sorted_frequency_pairs(k, 2) >= 1190.865 && sorted_frequency_pairs(k, 2) <= 1227.135
        detected_numbers(k) = 7;
        
    elseif  sorted_frequency_pairs(k, 1) >= 839.22 && sorted_frequency_pairs(k, 1) <= 864.78  && sorted_frequency_pairs(k, 2) >= 1315.96 && sorted_frequency_pairs(k, 2) <= 1356.04
        detected_numbers(k) = 8;
        
    elseif sorted_frequency_pairs(k, 1) >= 839.22 && sorted_frequency_pairs(k, 1) <= 864.78  && sorted_frequency_pairs(k, 2) >= 1454.845 && sorted_frequency_pairs(k, 2) <= 1499.155
        detected_numbers(k) = 9;
        
    elseif  sorted_frequency_pairs(k, 1) >= 839.22 && sorted_frequency_pairs(k, 1) <= 864.78  && sorted_frequency_pairs(k, 2) >= 1612.445 && sorted_frequency_pairs(k, 2) <= 1661.555
        detected_numbers(k) = 'C';
        
    elseif sorted_frequency_pairs(k, 1) >= 926.885 && sorted_frequency_pairs(k, 1) <= 955.115 && sorted_frequency_pairs(k, 2) >= 1190.865 && sorted_frequency_pairs(k, 2) <= 1227.135
        detected_numbers(k) = '*';
        
    elseif  sorted_frequency_pairs(k, 1) >= 926.885 && sorted_frequency_pairs(k, 1) <= 955.115  && sorted_frequency_pairs(k, 2) >= 1315.96 && sorted_frequency_pairs(k, 2) <= 1356.04
        detected_numbers(k) = 0;
        
    elseif sorted_frequency_pairs(k, 1) >= 926.885 && sorted_frequency_pairs(k, 1) <= 955.115  && sorted_frequency_pairs(k, 2) >= 1454.845 && sorted_frequency_pairs(k, 2) <= 1499.155
        detected_numbers(k)= '#';
        
    elseif  sorted_frequency_pairs(k, 1) >= 926.885  && sorted_frequency_pairs(k, 1) <= 955.115 && sorted_frequency_pairs(k, 2) >= 1612.445 && sorted_frequency_pairs(k, 2) <= 1661.555
        detected_numbers(k) = 'D';
        
    end
end
% Play the sound signal
sound(audio_signal)
% Display the detected numbers in the order they were dialed/processed.
detected_numbers'
