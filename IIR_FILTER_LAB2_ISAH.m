clc; % Clears the Command Window


% Read audio file and extract signal data and sampling frequency
[originalSignal, samplingFreq] = audioread("Don_Giovanni_1.wav");

% Calculate sampling period
samplingPeriod = 1 / samplingFreq;

% Determine signal length
signalLength = length(originalSignal);

% Create time vector corresponding to the signal
timeVector = 0:samplingPeriod:(signalLength-1)*samplingPeriod;

% Play the initial audio signal
sound(originalSignal, samplingFreq);

% Calculate the duration of the original signal for playback delay
originalDuration = length(originalSignal) / samplingFreq;

% Pause to allow the original signal to play before proceeding
pause(originalDuration + 1); % Add extra time to ensure playback completion


% Plot the first 100 samples of the signal against the time vector
figure;
plot(timeVector(1:1000), originalSignal(1:1000));
xlabel('Time (s)');
ylabel('Amplitude');

%% Spectral analysis: Compute frequency vector and perform FFT
frequencyStep = samplingFreq / signalLength;
frequencyVector = 0:frequencyStep:samplingFreq - frequencyStep;
signalFFT = abs(fft(originalSignal));

%% Plot the magnitude spectrum against frequency
figure;
plot(frequencyVector, signalFFT);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum');
grid on;

%% Perform FFT shift for centered frequency spectrum
frequencyVectorShift = -samplingFreq/2:frequencyStep:samplingFreq/2 - frequencyStep;
signalFFTShifted = abs(fftshift(fft(originalSignal)));

%% Plot the shifted magnitude spectrum against shifted frequency
figure;
plot(frequencyVectorShift, signalFFTShifted);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum with FFT Shift');
grid on;

%% Locate the first peak frequency
[maxPeak1, position1] = max(signalFFT(1:signalLength/2));
peakFrequency1 = frequencyVector(position1);

%% Calculate filter coefficients for the first peak frequency
zero1 = exp(1i*2*pi*peakFrequency1/samplingFreq);
zeroConj1 = conj(zero1);
pole1 = 0.99 * zero1;
poleConj1 = conj(pole1);

%% Define numerator and denominator coefficients for the filter transfer function
b0 = 1;
b1 = -(zero1 + zeroConj1);
b2 = zero1 * zeroConj1;
a0 = 1;
a1 = -(pole1 + poleConj1);
a2 = pole1 * poleConj1;

%% Create Z-plane plot for the first filter
figure;
zplane([zero1; zeroConj1], [pole1; poleConj1]);
title('Zeros and Poles Placement (Filter 1)');
legend('Zeros', 'Poles');
grid on;


%% Apply the first filter in the temporal domain
filteredSignal1 = filter([b0 b1 b2], [a0 a1 a2], originalSignal);


%% Compute FFT of the filtered signal 1
filteredSignalFFT1 = abs(fft(filteredSignal1));


%% Plot the magnitude spectrum of filtered signal 1
figure;
plot(frequencyVector, filteredSignalFFT1);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Filtered Signal 1');
axis([0 samplingFreq 0 100000]);
grid on;


%% Locate the second peak frequency in the filtered spectrum
[maxPeak2, position2] = max(filteredSignalFFT1(1:signalLength/2));
peakFrequency2 = frequencyVector(position2);


%% Calculate filter coefficients for the second peak frequency
zero2 = exp(1i*2*pi*peakFrequency2/samplingFreq);
zeroConj2 = conj(zero2);
pole2 = 0.99 * zero2;
poleConj2 = conj(pole2);


%% Define numerator and denominator coefficients for the second filter transfer function
b0_2 = 1;
b1_2 = -(zero2 + zeroConj2);
b2_2 = zero2 * zeroConj2;
a0_2 = 1;
a1_2 = -(pole2 + poleConj2);
a2_2 = pole2 * poleConj2;


%% Create Z-plane plot for the second filter
figure;
zplane([zero2; zeroConj2], [pole2; poleConj2]);
title('Zeros and Poles Placement (Filter 2)');
legend('Zeros', 'Poles');
grid on;


%% Apply the second filter in the temporal domain
filteredSignal2 = filter([b0_2 b1_2 b2_2], [a0_2 a1_2 a2_2], filteredSignal1);


%% Compute FFT of the final filtered signal 2
filteredSignalFFT2 = abs(fft(filteredSignal2));


%% Plot the magnitude spectrum of the final filtered signal 2
figure;
plot(frequencyVector, filteredSignalFFT2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Final Filtered Signal');
grid on;


%% Play the final filtered audio signal
sound(filteredSignal2, samplingFreq);
