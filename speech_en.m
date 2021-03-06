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
SNR = (abs(sum(Pyy./Pnn-1,1)));

figure(1)
plot(10*log10(SNR));
title('SNR')
axis([0 size(Y,2) -20 40])

S_ps = Spectral_Subtraction(Pyy, Pnn, Y);
S_w = Wiener1(Pyy, Pnn, Y);

s_out1 = stift(S_ps, win, l, o, 1, fs);
s_out2 = stift(S_w, win, l, o, 1, fs);
%snr_op = snr(s_out, fs)

%sound(s_out2, fs)
figure(2);
subplot(4,1,1)
plot(s(ind));
axis([0 70000 -0.5 0.5])
title('Clean Speech')
subplot(4,1,2)
plot(y);
axis([0 70000 -0.5 0.5])
title('Noisy Speech')
subplot(4,1,3)
plot(s_out1);
axis([0 70000 -0.5 0.5])
title('Estimated Clean Speech 1')
subplot(4,1,4)
plot(s_out2);
axis([0 70000 -0.5 0.5])
title('Estimated Clean Speech 2')