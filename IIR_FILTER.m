clc %% clears the Command Window, which is the window. 
close all %% close all open figure windows.


[signal, samp_freq] = audioread("Don_Giovanni_1.wav"); %%read the audio file (Don_Giovanni_1.wav)and store its waveform data in the "signal" variable and its sampling frequency in the "samp_freq" variable.

samp_period = 1/samp_freq; %%calculates the sampling period (time between consecutive samples) based on the sampling frequency.

len = length(signal); %%returns the number of elements (signal smaples) in the input signal and stores the value in len.

time_vec = 0:samp_period:(len-1)*samp_period; %%creates a time vector that corresponds to a given signal with a specified sampling period with 0 as the first value and (len-1)*samp_period at an interval of samp_period


figure %%creates a new figure window for displaying plots or graphics.
plot(time_vec(1:100), signal(1:100))%%plots a signal against its corresponding time vector, showing only the first 100 samples of the signal.
xlabel('time vector') %%labels the x-axis "time vector".
ylabel('Amplitude') %%labels the y-axis "Amplitude".

%% Spectral Analysis 
% Start by defining frequency vector
freq_step =samp_freq/len; % Spread the total number of samples equally from 0 to Sampling Frequency to determine the frequency resolution of the signal, which is the smallest change in frequency that can be detected

freq_vec = 0:freq_step:samp_freq-freq_step; %%generates a vector of frequencies that corresponds to a discrete Fourier transform (DFT) of a signal

SIGNAL = abs(fft(signal));%% computes the absolute value of the Discrete Fourier Transform (DFT) of a signal in MATLAB and stores the result in a variable called SIGNAL

figure %%creates a new figure window for displaying plots or graphics.
plot(freq_vec, SIGNAL) %%plots the freq_vec against the SIGNAL corresponding time vector.
xlabel('Frequency') %%labels the x-axis "Frequency".
ylabel('Magnitude') %%labels the y-axis "Magnitude".
grid on %%adds a grid to the plot.


%% Spectral Analysis using FFTSHIFT
freq_vec_shift = -samp_freq/2:freq_step:samp_freq/2-freq_step;
SIGNAL_shift = abs(fftshift(fft(signal))); 

%Figure
plot(freq_vec_shift, SIGNAL_shift)  %%plots the freq_vec_shift against the SIGNAL_shift corresponding time vector.
xlabel('Frequency') %%labels the x-axis "Frequency".
ylabel('Magnitude') %%labels the y-axis "Magnitude".
title('FFT SHIFT') %%adds a title to a MATLAB plot. The title of the plot is "FFT SHIFT"
grid on %%%%adds a grid to the plot.


%% LOCATING THE FIRST PEAK FREQUENCY
% The sampling freq should be greater than twice the max signal freq
[peak_1, position_1] = max(SIGNAL(1:len/2)); %%computes the maximum value and its position in a vector SIGNAL that represents the Discrete Fourier Transform (DFT) of a signal

freq_peak_1 = freq_vec(position_1-1); %%computes the frequency corresponding to the largest amplitude frequency component in a signal. 


%%
% filter zeros and poles for filter 1
zero_1 = exp(1i*2*pi*freq_peak_1/samp_freq);
zero_conj_1 = conj(zero_1);

%%
pole_1 = 0.99*zero_1;
pole_conj_1 = conj(pole_1);

%% Co-efficients of the numerator of filter transfer function
b0=1;
b1=-(zero_1+zero_conj_1);
b2=zero_1*zero_conj_1;

%%  Co-efficients of the denominator of filter transfer function
a0=1;
a1=-(pole_1+pole_conj_1);
a2= pole_1*pole_conj_1;

B1=[b0 b1 b2];
A1=[a0 a1 a2];

%% Z-Plane for the first filter
Z1=[zero_1;zero_conj_1];
P1=[pole_1;pole_conj_1];

figure
zplane(Z1,P1)
title ('Zeros and poles placment')
legend('zeros','poles')
grid on

%% Filtering is done in temporal domain
signal_filtered_1 = filter(B1,A1,signal);

SIGNAL_FILTERED_1 = abs(fft(signal_filtered_1));

figure
plot(freq_vec, SIGNAL_FILTERED_1)
xlabel('Frequency')
ylabel('Magnitude')
axis([0 samp_freq 0 100000]) %%sets the axis limits for the current plot.
grid on

%% LOCATING THE SECOND PEAK FREQUENCY
% The sampling freq should be greater than twice the max signal freq
[peak_2, position_2] = max(SIGNAL_FILTERED_1(1:len/2));

freq_peak_2 = freq_vec(position_2-1);

%% filter zeros and poles for filter 2
zero_2 = exp(1i*2*pi*freq_peak_2/samp_freq);
zero_conj_2 = conj(zero_2);

%%
pole_2 = 0.99*zero_2;
pole_conj_2 = conj(pole_2);

%% Co-efficients of the numerator of filter transfer function
b0_2=1;
b1_2=-(zero_2+zero_conj_2);
b2_2=zero_2*zero_conj_2;

%%  Co-efficients of the denominator of filter transfer function
a0_2=1;
a1_2=-(pole_2+pole_conj_2);
a2_2= pole_2*pole_conj_2;

B2=[b0_2 b1_2 b2_2];
A2=[a0_2 a1_2 a2_2];

%% Z-Plane for the first filter
Z2=[zero_2;zero_conj_2];
P2=[pole_2;pole_conj_2];

figure
zplane(Z2,P2)
title ('Zeros and poles placment')
legend('zeros','poles')
grid on

%% Filtering is done in temporal domain
signal_filtered_2 = filter(B2,A2,signal_filtered_1);

SIGNAL_FILTERED_2 = abs(fft(signal_filtered_2));

figure
plot(freq_vec, SIGNAL_FILTERED_2)
xlabel('Frequency')
ylabel('Magnitude')
grid on
%%
sound(signal_filtered_2, samp_freq)