clear all
close all

s = audioread('AudioFiles/clean_speech.wav');
type = 4;
db = 20;
n = genNoise(type, db, length(s));

ind = 1:length(s)/6;
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
[n_mmse,c_mmse] = mmse_noise_tracker(Pyy,Y);
n_mmse_spp = mmse_noise_tracker_spp(Pyy);
[~,~,n_vad] = vad_noise_tracker(y);

%% Variance

true_varw = zeros(size(Pnn,2),1);
true_varw(1) = mean(Pnn(:,1));
alpha = 0.8;
for i = 2:size(Pnn,2)
    true_varw(i) = mean(Pnn(:,i));
end

figure
t = linspace(0,size(Y,2), size(Y,2));%length(s)/fs);
plot(t,10*log10(true_varw))
hold on
plot(t,10*log10(mean(n_mmse_spp)),'-.')
hold on
plot(t,10*log10(mean(n_mmse)),'-')
hold on
plot(t,10*log10(mean(n_vad)),':')
xlabel('frame')
ylabel('\sigma^2_w (db)')
legend('True Noise Variance','MMSE-SPP Noise Variance','MMSE Noise Variance','VAD Noise Variance')
%axis([0 size(Y,2) -32 -17])
%% Clean Speech Estimation
c_vad = Wiener1(Pyy,n_vad,Y);
c_mmse_spp = Wiener1(Pyy,n_mmse_Y);

