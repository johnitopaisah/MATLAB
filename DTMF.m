% Read the audio file
[y, Fs] = audioread("0123456789.wav");
sound(y);

% Compute the sampling period Ts
Ts = 1/Fs;

% Compute the signal length of the number of samples in the signal
N = length(y);

% Create the time vector for the audio signal
t = 0:Ts:(N-1)*Ts;

% Plot the original audio signal in the time domain
figure(1);
plot(t, y);
xlabel('Time in seconds');
ylabel('Sampled Signal');
title('Given Signal in Time Domai');

% The frequency range of the audio signal using the Nyquist frequency and
% the number of samples
f = Fs/2;
M = 256; % 256 is choosen for an optimum frequency resolution forr the 
% spectrum

% Create a linearly spaced freqency victor from 0hz to sampling frequency
f1 = linspace(0, Fs, M+1);
% Remove the last element from the frequency vector and discard the Nyquist
% frequency point from the vector.
f1 = f1(1:end-1);
%  Convert the frequencies from the Nyquist frequency to the normalized 
% range of 0 to 1, where 1 represents the Nyquist frequency.
f2 = f1/Fs;
% Shift the frequencies in ht e vector f1 to be centered around zero
f1 = f1 - f;

% Compute the FFT of the signal with the resolution M to examine the
% magnitude and phase information in the resulting spectrum
Y = fft(y, M);

figure (2);
plot(f2, abs(Y)) % Plot the magnitude spectrum of the original audio signal
% in frequency domain
title('Spectrum of the original signal before filtering');
xlabel('Frequency of the Signal');
ylabel('Spectrum of the Signal');
grid on;

% Spectrum analysis of each window in 65ms
D = 65e-3;
% Determine the number of samples for each window
nn = round(D/Ts);

% Reshaping the signal in intervals of the window size
fint = 1/N; % interval between frequency samples
% Creating frequency vector containing th frequencies corresponding to each
% frequency
f = 0:fint:(N-1)*fint;
% Sample point in the Fourier Transform

% Reshaping the signal in intervals of the window size
nb = ceil(N/nn)*nn - N; % Determine the number of zeros that need to be
% added to the signal 
ye1 = [y.' zeros(1, nb)]; % Concatenate the original signal y with a vector 
% of zeros of length nb
y_tab = reshape(ye1', nn, [])'; % Reshape the signal by transposing ye1
y_tab2 = zeros(size(y_tab)); % Initialize a matrix of zeros with the same
% size as y_tab

y_tab1 = [t zeros(1, nb)]; % Concatenate the time victor t with 
% a vector of zeros of lenth nb
y_tab2 = reshape(y_tab1', nn, [])'; % Restore the original orientation and 
% has the same dimensions as the reshaped signal matrix y_tab

% Minimum power initialiization
power_min = 20;

% Initialization of the victor
Tone = zeros(size(y_tab, 1), 1); % Create a placeholder for storing results
% related to the processing of the signal sigments

% Iterate over the rows of the y_row matrix
for k = 1:size(y_tab, 1)

    power = sum(y_tab(k, :).^2)/nn; % Power of the signal
    power_dB = 10*log10(power/10^-3); % Power of the signal in dB

    Tone(k) = power_dB > power_min; % Select all windows with power 
    % greater than 20 dB
end

% Compute the numberr of windows or signal sigments that have power levels
% are greater than the threshold of 20dB
s1 = sum(Tone);

power = zeros(size(y_tab, 1), 1); % Placeholder used to shtore the 
% calculated power of each signal segment
power_dB = zeros(size(y_tab, 1), 1); % Used to store the power values 
% converted to decibels (dB)
Tone2 = zeros(s1, 2); % Used to store the frequencies obtained from the 
% signal segments that have power levels above the specified threshold

% Parameters for cheby1() filter function
aa = 7;
bb = 10;
cc = 0.55;
dd = 0.93;

% Iterating over rows of y_tab
for k = 1:size(y_tab, 1)
    if k == 1 
        % Calculate the Fast Fourier Transform (FFT) of the k-th row of
        % y_tab using N points.
        Y_TAB1 = fft(y_tab(k, :), N);
        % Compute the magnitude of the FFT retuslt
        Y_TAB1_mag = abs(Y_TAB1);
        %  Find the peaks in the magnitude spectrum of the positive
        %  frequencies
        [peak1, loc1] = findpeaks(Y_TAB1_mag(1:(N+1)/2), 'SortStr','descend');
        % Assign the frequency corresponding to the highest peak to the
        % first column of Tone2
        Tone2(k,1) = f(loc1(1)); % to get the location of the first frequency
        N1 = Tone2(k,1); % Stores the frequency value

        % Design a Chebyshev Type I filter with the specified paramters
        % based on the frequency range [cc*N1 dd*N1]
        [b,a] = cheby1(aa, bb, [cc*cc*N1 dd*N1], 'stop');
        % Apply the designed filter b and a to the k-th row of y_tab and
        % stores the filtered result
        y_tab2(k,:) = filter(b,a,y_tab(k,:));
        % Calculate the FFT of the filterred signal
        Y_TAB2 = fft(y_tab2(k,1),N);
        % Calculate the magnitude spectrum of the filtered signal.
        Y_TAB2_mag = abs(Y_TAB2);
        % Finds the peaks in the magniitude spectrum of the positive
        % frequencies
        [peak2, loc2] = findpeaks(Y_TAB2_mag(1:(N+1)/2), 'SortStr','descend');
        % Assign the frequency corresponding to the highest peak to the
        % second colunm of Tone2
        Tone2(k, 2) = f(loc2(1));
    else
        %  For iterations other than the first one
        power(k-1) = sum(y_tab(k-1,:).^2)/nn; % Power of the signal
        power_dB(k-1) = (10*log20(power(k-1)/10.^-3)); % Power of the signal in dB

        power(k) = sum(y_tab(k-1,:).^2)/nn; % Power of the signal
        power_dB(k) = (10*log20(power(k-1)/10.^-3)); % Power of the signal in dB

        if power_dB(k-1) < power_min && power_dB(k) > power_min % To avoid taking duplicate freqencies
            
            Y_TAB1 = fft(y_tab(k,:), N); % Calculate the FFT of the k-th row of y_tab using N points.

            Y_TAB1_mag = abs(Y_TAB1);
            % Find the peaks in the magnitude spectrum of the positive
            % frequencies up to Nyquist frequency.
            [peak1, loc1] = findpeaks(Y_TAB1_mag(1:(N+1)/2), 'SortStr','descend');
            % Assign the frequency corresponding to the highest peak to the
            % first column of Tone2.
            Tone2(k,1) = f(loc1(1));

            N1 = Tone2(k,1); % Stores the frequency value in N1.
            [b,a] = cheby1(aa,bb,[cc*N1 dd*N1], 'stop'); % To design the filter

            y_tab2(k,:) = filter(b,a, y_tab(k,:)); % To apply the filter
            Y_TAB2 = fft(y_tab2(k,:), N); % Compute the FFT of the filtered signal
            Y_TAB2_mag = abs(Y_TAB2); % Calculate the magnitude spectrum of the filtered signal
            % Find the peaks in the magnitude spectrum of the positive
            % frequencies up to Nyquist freqency
            [peak2, loc2] = findpeaks(Y_TAB2_mag(1:(N+1)/2), 'SortStr','descend');
            Tone2(k,2) = f(loc2(1)); % To get the location of the second frequency

        end
    end
end

% The first and second column were sorted to make the first column to be 
% the low frequency and second column the high frequencies
Tone3 = reshape(nonzeros(Tone2), [], 2)*Fs; % Obtain the frequencies in Hz
Tone4 = sort(Tone3.','ascend'); % Transpose Tone3 matrix and sorts it in 
% ascending order resulting
Tone5 = Tone4'; % Restore the original arrangement of frequency pairs
number = zeros(length(Tone5), 1); % To store the calculated values (numbers dailed) in the subsequent code

% To compare the frequencies in the signal to the DTMF, considering the
% 1.5% frequency tolerance

for k = 1:length(Tone5)

    if ~isempty(Tone5) && size(Tone5, 2) >= 2 && Tone5(k,1) >= 686.545 && Tone5(k,1) <= 707.455 && Tone5(k,2) >= 1990.865 && Tone5(k,2) <= 1227.135
        number(k) = 1;
    elseif Tone5(k,1) >= 686.545 && Tone5(k,1) <= 707.545 && Tone5(k,2) >= 1315.96 && Tone5(k,2) <= 1356.04
        number(k) = 2;
    elseif Tone5(k,1) >= 686.545 && Tone5(k,1) <= 707.545 && Tone5(k,2) >= 1415.445 && Tone5(k,2) <= 1499.155
        number(k) = 3;
    elseif Tone5(k,1) >= 686.545 && Tone5(k,1) <= 707.545 && Tone5(k,2) >= 1612.445 && Tone5(k,2) <= 1661.555
        number(k) = 'A';

    elseif Tone5(k,1) >= 758.45 && Tone5(k,1) <= 781.55 && Tone5(k,2) >= 1190.865 && Tone5(k,2) <= 1227.135
        number(k) = 4;
    elseif Tone5(k,1) >= 758.45 && Tone5(k,1) <= 781.55 && Tone5(k,2) >= 1315.96 && Tone5(k,2) <= 1356.04
        number(k) = 5;
    elseif Tone5(k,1) >= 758.45 && Tone5(k,1) <= 781.55 && Tone5(k,2) >= 1454.845 && Tone5(k,2) <= 1499.155
        number(k) = 6;
    elseif Tone5(k,1) >= 758.45 && Tone5(k,1) <= 781.55 && Tone5(k,2) >= 1612.445 && Tone5(k,2) <= 1661.555
        number(k) = 'B';

    elseif Tone5(k,1) >= 839.22 && Tone5(k,1) <= 864.78 && Tone5(k,2) >= 1190.865 && Tone5(k,2) <= 1227.135
        number(k) = 7;
    elseif Tone5(k,1) >= 839.22 && Tone5(k,1) <= 864.78 && Tone5(k,2) >= 1315.96 && Tone5(k,2) <= 1356.04
        number(k) = 8;
    elseif Tone5(k,1) >= 839.22 && Tone5(k,1) <= 864.78 && Tone5(k,2) >= 1454.845 && Tone5(k,2) <= 1499.155
        number(k) = 9;
    elseif Tone5(k,1) >= 839.22 && Tone5(k,1) <= 864.78 && Tone5(k,2) >= 1612.445 && Tone5(k,2) <= 1661.555
        number(k) = 'C';

    elseif Tone5(k,1) >= 926.885 && Tone5(k,1) <= 955.115 && Tone5(k,2) >= 1190.865 && Tone5(k,2) <= 1227.135
        number(k) = '*';
    elseif Tone5(k,1) >= 926.885 && Tone5(k,1) <= 955.115 && Tone5(k,2) >= 1315.96 && Tone5(k,2) <= 1356.04
        number(k) = 0;
    elseif Tone5(k,1) >= 926.885 && Tone5(k,1) <= 955.115 && Tone5(k,2) >= 1354.854 && Tone5(k,2) <= 1499.155
        number(k) = '#';
    elseif Tone5(k,1) >= 926.885 && Tone5(k,1) <= 955.115 && Tone5(k,2) >= 1612.445 && Tone5(k,2) <= 1661.555
        number(k) = 'D';
    end
end
