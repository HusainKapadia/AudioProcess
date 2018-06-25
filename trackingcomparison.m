clear all
close all

s = audioread('AudioFiles/clean_speech.wav');
type = 2;
db = 20;
n = genNoise(type, db, length(s));

ind = 1:96240;
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
[n_mmse,c_mmse_fft] = mmse_noise_tracker(Pyy,Y,'Wiener');
[n_mmse_spp,c_mmse_spp_fft] = mmse_noise_tracker_spp(Pyy,Y,'Wiener');
[n_vad,c_vad_fft] = vad_noise_tracker(Pyy,Y,'Wiener');

%% Variance

true_varw = zeros(size(Pyy,2),1);
true_varw(1) = mean(Pyy(:,1));
alpha = 0.8;
for i = 2:size(Pyy,2)
    true_varw(i) = mean(Pyy(:,i));
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

err_mmse = mean(mean(abs(10*log(Pnn./n_mmse))));
err_mmse_spp = mean(mean(abs(10*log(Pnn./n_mmse_spp))));
err_vad = mean(mean(abs(10*log(Pnn./n_vad))));
%% Clean Speech Estimation
%c_vad_fft = Wiener1(Pyy,n_vad,Y);
c_vad = stift(c_vad_fft, win, l, o, 1, fs);

%c_mmse_spp_fft = Wiener1(Pyy,n_mmse_spp,Y);
c_mmse_spp = stift(c_mmse_spp_fft, win, l, o, 1, fs);

%c_mmse_fft = Wiener1(Pyy,n_mmse,Y);
c_mmse=stift(c_mmse_fft, win, l, o, 1, fs);
%% SSNR
ssnr_vad = ssnr(s(ind),c_vad,fs);
ssnr_mmse_spp = ssnr(s(ind),c_mmse_spp,fs);
ssnr_mmse = ssnr(s(ind),c_mmse,fs);
ssnr_noisy = ssnr(s(ind),y(ind),fs);
%% STOI
stoi_vad = stoi(s(ind),c_vad,fs);
stoi_mmse_spp =stoi(s(ind),c_mmse_spp,fs);
stoi_mmse = stoi(s(ind),c_mmse,fs);
stoi_noisy = stoi(y(ind),s(ind),fs);