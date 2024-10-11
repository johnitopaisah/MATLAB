[y, Fs] = audioread("Don_Giovanni_1.wav");
N = length(y);
player = audioplayer(y, Fs);

t = fft(y, N);
Y = fftshift(t);
f = -Fs/2:Fs/N:(Fs/2-Fs/N);
Te = 1/Fs;

subplot(2, 1, 1);
plot(y);
title('Time domain waveform of the original audio signal');
xlabel('Time');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
plot(f, abs(t));
title('Original signal in Frequency-domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

[val, i] = max(Y);
