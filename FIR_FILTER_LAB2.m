clc;
%%clear all;
close all;

% Load audio signal
[audioSignal, samplingFreq] = audioread('Don_Giovanni_2.wav');
signalLength = length(audioSignal);
samplingPeriod = 1 / samplingFreq;
timeVector = 0:samplingPeriod:(signalLength-1)*samplingPeriod;
frequencyMax = samplingFreq / 2;
frequencyInterval = samplingFreq / signalLength;
frequencyVector = 0:frequencyInterval:(signalLength-1)*frequencyInterval;

% Compute FFT of the input audio signal
fftResult = fft(audioSignal);
fftMagnitude = abs(fftResult);

% Define parameters for the moving average filter
windowSize = 31;
windowFunction = gausswin(2*windowSize+1) ./ windowSize;

% Apply the digital filter to the input signal
filteredSignal = filter(windowFunction, 1, audioSignal);
sound(filteredSignal, samplingFreq);

% Compute FFT of the filtered signal
filteredFFTResult = fft(filteredSignal);
filteredFFTMagnitude = abs(filteredFFTResult);

% Define parameters for sinc functions
sincArg1 = pi * samplingFreq * 3 * samplingPeriod * timeVector;
sincArg2 = pi * samplingFreq * samplingPeriod * timeVector;
sincResult1 = sinc(sincArg1);
sincResult2 = sinc(sincArg2);

% Compute G by dividing sincResult1 with sincResult2
G = sincResult1 ./ sincResult2;

% Compute FFT of G
GFFTResult = fft(G);
GFFTMagnitude = abs(GFFTResult);
phaseAngles = -2 * pi * frequencyVector * (signalLength-1) / 2 * samplingPeriod;

% Plot time domain analysis
figure(1);
plot(timeVector, audioSignal, 'r');
grid on;
hold on;
plot(timeVector, filteredSignal, 'b');
title('Time Domain Analysis N=31');
xlabel('Time in seconds');
ylabel('Amplitude');
legend('Before Filtering Time', 'After Filtering Time');

% Plot frequency domain analysis
figure(2);
semilogy(frequencyVector, fftMagnitude, 'y');
grid on;
hold on;
semilogy(frequencyVector, filteredFFTMagnitude, 'g');
title('Spectrum of Signal N=31');
xlabel('Frequency in Hz');
ylabel('Amplitude in log scale');
legend('Before Filtering Freq', 'After Filtering Freq');
