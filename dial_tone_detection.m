% Dial Tone Detection
%% clear the command window, remove all variables from the workspace, and close all open figures
clc;
clear;
close all;

% Read the audio signal
[y, Fs] = audioread('numero2.wav');

% Preprocessing
y = resample(y, 8000, Fs); % Resample to 8kHz
y = y / max(abs(y)); % Normalize

% Dial Tone Detection Parameters
dialToneFrequency = 680; % Dial tone frequency in Hz
bandwidth = 10; % Bandwidth for the bandpass filter in Hz
threshold = 0.1; % Threshold for dial tone detection

% Apply bandpass filter
%fNorm = dialToneFrequency / (8000 / 2);

%new
% Convert frequency and bandwidth to normalized values
fNorm = dialToneFrequency / (Fs / 2);
bwNorm = bandwidth / (Fs / 2);


% Apply bandpass filter
[b, a] = butter(2, [fNorm - bwNorm/2, fNorm + bwNorm/2], 'bandpass');
filteredSignal = filter(b, a, y);

% Envelope detection
envelope = abs(filteredSignal);
%smoothedEnvelope = smooth(envelope, 200); % Apply moving average or low-pass filter

% Apply moving average filter
windowSize = 200;
movingAvgFilter = ones(1, windowSize) / windowSize;
smoothedEnvelope = conv(envelope, movingAvgFilter, 'same');

% Thresholding
dialToneDetected = smoothedEnvelope > threshold;

% DTMF Decoding

% Define DTMF frequencies
frequencies = [697 770 852 941; 1209 1336 1477 1633];

% Define DTMF digits
digits = ['1' '2' '3' 'A'; '4' '5' '6' 'B'; '7' '8' '9' 'C'; '*' '0' '#' 'D'];

% DTMF Decoding Parameters
digitDuration = 0.065; % Duration of a dialed digit in seconds
digitGap = 0.08; % Delay between two dialed digits in seconds

% Perform DTMF decoding
samplesPerDigit = digitDuration * 8000;
samplesPerGap = digitGap * 8000;

digitStart = 1;
phoneNumber = '';

while digitStart <= length(y)
    % Check if dial tone detected
    if dialToneDetected(digitStart)
        % Extract digit samples
        digitSamples = y(digitStart : digitStart + samplesPerDigit - 1);
        
        % Apply Goertzel algorithm for frequency detection
        for i = 1:size(frequencies, 1)
            f1 = frequencies(i, 1);
            f2 = frequencies(:, 2:end);
            power = goertzel(digitSamples, f1) + goertzel(digitSamples, f2);
            
            % Check if power exceeds threshold
            if power > threshold
                phoneNumber = [phoneNumber digits(i)];
                break;
            end
        end
        
        digitStart = digitStart + samplesPerDigit;
    else
        digitStart = digitStart + samplesPerGap;
    end
end

disp(['Decoded Phone Number: ' phoneNumber]);
%sound(y)
