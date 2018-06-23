clc;
close all;

s = audioread('AudioFiles/clean_speech.wav');
type = 2;
db = 20;
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
varw_hat(:,1) = Pyy(:,1);
smoothed_varw_hat = zeros(size(Pyy));
%lrt = zeros(1, size(Pyy,2));

%Update noise using safety net
time_interval = 800; %samples
time_frame = floor(time_interval/l); %safety frame
alpha = 0.99;
beta = 0.9;
S_hat = 0;
[lrt,silent_frames,varw_hat] = vad_noise_tracker(y);
t = linspace(0,length(y)/fs, length(ind));%length(s)/fs);
figure
plot(t,s(ind));
hold on
plot(t,silent_frames,'linewidth',6)
xlabel('time (s)')
ylabel('Amplitutude')
legend('Clean Speech','Silent Periods')
figure
plot(t,s(ind));
hold on
plot(t,lrt)
xlabel('time (s)')
ylabel('log\lambda')
legend('Clean Speech','log\lambda')
true_varw = zeros(size(Pnn,2),1);
true_varw(1) = mean(Pnn(:,1));
for i = 2:size(Pnn,2)
    true_varw(i) = beta*mean(Pnn(:,i-1))+(1-beta)*mean(Pnn(:,i));
end

figure
t = linspace(0,size(Y,2), size(Y,2));%length(s)/fs);
plot(t,10*log10(true_varw))
hold on
plot(t,10*log10(mean(varw_hat)))
xlabel('frame')
ylabel('\sigma^2_w (db)')
legend('True Noise Variance','Vad Based Noise Variance')
%axis([0 size(Y,2) db-5 db+20])

% figure
% plot(L)
% xlabel('Number of frames')
% ylabel('Likelihood Ratio')

S_ps = Spectral_Subtraction(Pyy, varw_hat, Y);
s_out1 = stift(S_ps, win, l, o, 1, fs);

S_w = Wiener1(Pyy, varw_hat, Y);
s_out2 = stift(S_ps, win, l, o, 1, fs);

%sound(s_out1, fs)
%sound(s_out2, fs)
