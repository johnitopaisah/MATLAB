clc;
clear all;
close all;
%Synthesis of audio signal
[Y, Fe] = audioread('Don_Giovanni_2.wav');
N = length(Y);
Te = 1/Fe;
t=0:Te:(N-1)*Te;
f= Fe/2;
fint = Fe/N;
f = 0:fint:(N-1)*fint;
X=fft(Y);
X_mag = abs(X);
window = 31; %%represent the size of a moving average window or the length of a window function used for spectral analysis.
% The walue of the window can also be chnaged to 3 and 11 as indicated in the question.
h = gausswin(2*window+1)./window;%%The value of h is calculated using the gausswin function, which creates a Gaussian window of length 2*window+1
y = filter(h,1,Y); %%applying a digital filter to the input signal Y.
sound(y, Fe)
Y=fft(y);%%performs a Fast Fourier Transform (FFT) on the filtered signal y.
Y_mag = abs(Y);
a = sinc(pi*Fe*3*Te*t); %%The value of a is calculated using the sinc function, which generates a sinc function of the form sinc(x) = sin(pi*x)/(pi*x)
b = sinc(pi*Fe*Te*t);
g =a./b;
G = fft(g);%%performs a Fast Fourier Transform (FFT) on the signal g.
G_mag = abs(G);%%calculates the magnitude of the complex values in the variable G, and assigns the result to the variable G_mag
ang = -2*pi*f*(N-1)/2*Te;%%calculates a vector ang, which represents the phase angles of each frequency component in the FFT output G.
fig=1;
figure(fig) %%creates a new figure window and sets it as the current figure for plotting, using the value of the variable fig as the figure identifier.

plot(t, Y, 'r')%%creates a plot of the signal Y as a function of time t, with the color of the plot set to red.
grid
hold on %%retain the current plot while adding new plots to the same figure.
plot(t, y, 'b')%%creates a plot of the signal y as a function of time t, with the color of the plot set to blue.
title('TIme Domain Analysis N=31')
xlabel('TIme in secs')
ylabel('Amplitude')
legend('Before Filtering Time','After Filtering Time')%%adds a legend to a plot, with two labels: "Before Filtering Time" and "After Filtering Time".
fig=fig+1;
figure(fig)
semilogy(f, X_mag, 'y')%%plots the magnitude of a Fourier transform X on a semilogarithmic scale, using a yellow line.
grid
hold on
semilogy(f, Y_mag, 'g')
title('Spectrum of Signal N=31')
xlabel('Frequency in Hz')
ylabel('Amplitude in log scale')
legend('Before Filtering Freq', 'After Filtering Freq')