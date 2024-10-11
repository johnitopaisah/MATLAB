%     Audio Signal Processing in MATLAB
%     
%     STUDENT_1 NAME: John Itopa ISAH - 63362
%     STUDENT_2 NAME: Blesso Danny Jeba Stephenen - 63291
%
%
%   Introduction:
%    This MATLAB script demonstrates basic audio signal processing techniques,
%    including Fast Fourier Transform (FFT), frequency manipulation, and inverse FFT.
%    It reads an audio signal from a WAV file, performs frequency manipulation, and
%    plays both the original and processed audio signals.
%
%    Usage:
%    1. Prerequisites: MATLAB software installed on your system.
%    2. Download the Code: Download the MATLAB script (`SignalSwap.m`) to your working directory.
%    3. Download Audio File: Ensure you have an audio file named `canal.wav` in the same directory as the MATLAB script.
%       You can replace it with your own audio file if desired.
%    4. Run the Script: Open MATLAB, navigate to the directory containing the script, and run the script by executing
%       the command `SignalSwap.m` in the MATLAB command window.
%    5. Listening to the Audio: The script will play the original audio signal first, followed by the processed audio signal.
%       Listen to both signals to observe the difference.
%    6. Visualizing Frequency Content: The script will also generate plots of the FFT magnitude before and after processing,
%       allowing you to visualize the frequency content of the signals.
%
%    Code Explanation:
%    - The script begins by clearing the command window, closing all open figures, and reading the audio signal from the file 'canal.wav'.
%    - It then plays the original audio signal using the `sound` function and waits for the playback to finish.
%    - Next, the script calculates the FFT of the original signal and performs frequency manipulation to swap high and low frequencies.
%    - The processed signal is obtained by converting the manipulated FFT back to the time domain using the inverse FFT.
%    - The processed signal is played using the `sound` function, and plots of the FFT magnitude before and after processing
%      are generated for visualization.
%
%    Files:
%    - `audio_processing.m`: MATLAB script for audio signal processing.
%    - `canal.wav`: Sample audio file used for demonstration.
%
%    Dependencies:
%    - MATLAB Signal Processing Toolbox (for FFT and inverse FFT).
%
%    Note:
%    - Ensure that the audio file 'canal.wav' is present in the same directory as the script.
%    - Adjust the file name and path if using a different audio file.



% Clear the command window, remove all variables from the workspace, and
% close all open figures
clc;
close all;

% Read the audio signal from the file 'canal.wav' and store the signal data
% in originalSignal and the sampling frequency in samplingFrequency.
[originalSignal, samplingFrequency] = audioread('canal.wav');

% Play the original sound
sound(originalSignal, samplingFrequency);

% Wait for the original sound to finish playing before processing
pause(length(originalSignal)/samplingFrequency);

% Establish the signal length or number of samples
numberOfSamples = length(originalSignal);

% Perform a Fast Fourier Transform (FFT) of the signal and store the results 
% in fftResult.
fftResult = fft(originalSignal);

% Calculate the magnitude of the FFT of the signal and store it in fftMagnitude.
fftMagnitude = abs(fftResult);

% Determine the length of the FFT
fftLength = length(fftResult);

% Create a new complex vector newFftResult from the FFT results with positive and
% negative frequencies
newFftResult(1) = fftResult(1);
newFftResult(fftLength/2) = fftResult(fftLength/2);
newFftResult(2:1:fftLength/2-1) = fftResult(fftLength/2-1:-1:2);
newFftResult(fftLength/2+1:1:fftLength) = fftResult(fftLength:-1:fftLength/2+1);

% Convert the processed signal back to the time domain using the inverse FFT
processedSignal = ifft(newFftResult);

% Play the processed sound
sound(real(processedSignal), samplingFrequency);

% Create a time vector for the audio signal
timeVector = (0:numberOfSamples-1) / samplingFrequency;

% Create a frequency vector using the Nyquist frequency and the number of
% samples
frequencyInterval = samplingFrequency/numberOfSamples;
frequencyVector = 0:frequencyInterval:(numberOfSamples-1)*frequencyInterval;

% Calculate the magnitude of the newFftResult vector and store it in newFftMagnitude
newFftMagnitude = abs(newFftResult);

% Create a new figure and plot the original signal's FFT magnitude
originalFig = 1;
figure(originalFig);
plot(fftMagnitude);
grid
title('Signal Before Processing');
xlabel('Frequency');
ylabel('Amplitude');

% Increment the figure counter and create a new figure to plot the
% processed signal's FFT magnitude
processedFig = originalFig + 1;
figure(processedFig);
plot(newFftMagnitude);
grid
title('Signal After Processing')
xlabel('Frequency');
ylabel('Amplitude');
