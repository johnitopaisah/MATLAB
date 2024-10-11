[y, Fs] = audioread("Don_Giovanni_1.wav");
player = audioplayer(y, Fs);
play(player) 
% sound(y, Fs)

N = length(y);
Y = fft(y);

frequencies = linspace(0, Fs/2, N/2+1);
magnitudeSpectrum = abs(Y(1:N/2+1));

plot(frequencies, magnitudeSpectrum);
xlabel('Frequencies (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum of the Audio Signal');