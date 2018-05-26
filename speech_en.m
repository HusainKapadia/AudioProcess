clear all;
close all;

s = audioread('AudioFiles/clean_speech.wav');
n = 0.5*audioread('AudioFiles/babble_noise.wav');
ind = 1:70000;
fs = 16000;
l = 15;
o = 60;
win = 4;
snr_ip = snr(s(ind), n(ind))

y = s(ind)+ n(ind);
snr_ip1 = snr(y, fs)

Y = stft(y, win, l, o, 1, fs);
N = stft(n(ind), win, l, o, 1, fs);

Pyy = Bartlett_P(Y, 8);
Pnn = Bartlett_P(N, 8);

figure(1)
subplot(2,1,1)
plot(Pyy(:,118));
title('Noisy Speech Periodogram, Pyy')
subplot(2,1,2)
plot(Pnn(:,118));
title('Noise Periodogram, Pnn')

S = Spectral_Subtraction(Pyy, Pnn, Y);

s_out = stift(S, win, l, o, 1, fs);
snr_op = snr(s_out, fs)

sound(s_out, fs)
figure(2);
subplot(3,1,1)
plot(s(ind));
axis([0 70000 -0.5 0.5])
title('Clean Speech')
subplot(3,1,2)
plot(y);
axis([0 70000 -0.5 0.5])
title('Noisy Speech')
subplot(3,1,3)
plot(s_out);
axis([0 70000 -0.5 0.5])
title('Estimated Clean Speech')