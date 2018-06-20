clear all;
close all;

s = audioread('AudioFiles/clean_speech.wav');
%n = 0.5*wgn(length(s),1,-80);
n = 0.5*audioread('AudioFiles/babble_noise.wav');
ind = 1:70000;
fs = 16000;
l = 15;
o = 60;
win = 4;
% snr_ip = snr(s(ind), n(ind))

y = s(ind)+ n(ind);
% snr_ip1 = snr(y, fs)

Y = stft(y, win, l, o, 1, fs);
N = stft(n(ind), win, l, o, 1, fs);
M = 8;
Pyy = Bartlett_P(Y, M);
Pnn = Bartlett_P(N, M);

%Initialize the noise estimate by assuming first eight frames are noise
varw_hat = zeros(size(Pyy));
varw_hat(:,1) = mean(Pyy(:,1:5),2);
[varw_hat,S_w] = mmse(Y);
true_varw = zeros(size(Pnn,2),1);
for i = 1:size(Pnn,2)
    true_varw(i) =mean(Pnn(:,i));
end

figure
t = linspace(0,size(Y,2), size(Y,2));%length(s)/fs);
plot(10*log10(true_varw))
hold on
plot(10*log10(mean(varw_hat)))
xlabel('time (s)')
ylabel('\sigma^2_w (db)')
legend('true noise variance','mmse estimated noise variance')

S_ps = Spectral_Subtraction(Pyy, varw_hat, Y);
s_out1 = stift(S_ps, win, l, o, 1, fs);

S_w2 = Wiener1(Pyy, varw_hat, Y);
s_out2 = stift(S_w, win, l, o, 1, fs);
err = sum(sum((abs(S_w2(:,1:720)) - abs(S_w(:,1:720))))).^2
%sound(s_out1, fs)
sound(s_out2, fs)
figure
subplot(3,1,1)
plot(s_out2);
axis([0 70000 -0.7 0.7])
title('Estimated Clean Speech 1')
subplot(3,1,2)
plot(y);
axis([0 70000 -1 1])
title('Noisy Speech')
subplot(3,1,3)
plot(s);
axis([0 70000 -0.7 0.7])
title('Clean Speech')