clc;
close all;

s = audioread('AudioFiles/clean_speech.wav');
type = 1;
db = 28;
n1 = genNoise(type, 25, length(s));
n2 = genNoise(2, db, length(s));
n3 = genNoise(3, db, length(s));
n4 = genNoise(4, 40, length(s));

fs = 16000;                     %sampling frequency
ind = 1:length(s);                  %length of speech
m = 4;                          %number of mics
l = 15;                         %frame length in ms
o = 60;                         %percent overlap

%% Generate Noisy Speech
c = 340;                        %speed of sound (m/s)
a = 40;                         %target direction angle (degrees)
x = 0.1*(1:m-1);                %distance between mics (m) 
tau = x.*fs.*cos(pi*a/180)./c;  %Time delay (s)

s = repmat(s(ind), 1, m);                      %clean speech
n = [n1(ind) n2(ind) n3(ind) n4(ind)];         %noise
    
%%STFT with overlap
S = stft(s, 4, l, o, m, fs);
N = stft(n, 4, l, o, m, fs);
Y = zeros(size(S));

S_new = permute(S, [1 3 2]);
N_new = permute(N, [1 3 2]);
Y_new = permute(Y, [1 3 2]);

d = [ones(1, size(Y,1)); exp(-1i*2*pi*(1:size(Y,1))*tau(1)/size(Y,1)); exp(-1i*2*pi*(1:size(Y,1))*tau(2)/size(Y,1)); exp(-1i*2*pi*(1:size(Y,1))*tau(3)/size(Y,1))].';          %steering vector
for k=1:size(S,2)
    Y_new(:,:,k) = d.*S_new(:,:,k) + N_new(:,:,k);
end

Y = permute(Y_new, [1 3 2]);
y = stift(Y, 3, l, o, m, fs);

%% Synthesize clean speech
Cw = zeros(m);

P = permute(Y_new, [2 1 3]);
for j = 1:size(Y,2)
    Z = P(:,:,j);
    %%Noise Covariance
    Cw = (j*Cw + corr(Y_new(:,:,j)))/(j+1);
    for k =1:size(Y,1)
        %Delay & Sum Beamformer
        %W = (1/m)*d(k,:).';
        %%MVDR Beamformer
        W = (pinv(Cw)*d(k,:).')/(d(k,:)*pinv(Cw)*d(k,:)');
        S_e(k,j) = W'*Z(:,k);
    end
end   

%%Single channel Wiener
M = 8;
Pse = Bartlett_P(S_e, M);
%[noise_psd1, S_o1] = mmse_noise_tracker(Pse, S_e, 'SS');
%[noise_psd2, S_o2] = vad_noise_tracker(Pse, S_e, 'Wiener');
[noise_psd3, S_o3] = mmse_noise_tracker_spp(Pse, S_e, 'SS');
%[noise_psd, S_o] = mmse_noise_tracker_spp(Pse, Y, 'Wiener');

%%STIFT with overlap add
s_e = stift(S_e, 4, l, o, 1, fs);
%s_o1 = stift(S_o1, 4, l, o, 1, fs);
%s_o2 = stift(S_o2, 4, l, o, 1, fs);
s_o3 = stift(S_o3, 4, l, o, 1, fs);

%% Plots
t = (1:length(y))./fs;
figure()
subplot(3,1,1)
plot(t, y), title('Noisy Speech');
axis([0 length(t)/fs -0.7 0.7]);
subplot(3,1,2)
plot(t, s_e), title('MVDR output');
axis([0 length(t)/fs -0.7 0.7]);
subplot(3,1,3)
plot(t, s_o3), title('Estimated Clean Speech');
axis([0 length(t)/fs -0.7 0.7]);
%% Evaluation
%SSNR1 = ssnr(s(1:length(s_o1),1), s_o1, fs)
%STOI1 = stoi(s(1:length(s_o1),1), s_o1, fs)
% 
% SSNR2 = ssnr(s(1:length(s_o2),1), s_o2, fs)
% STOI2 = stoi(s(1:length(s_o2),1), s_o2, fs)
% 
SSNR = ssnr(s(1:length(y),1), y(:,1), fs);
SSNR3 = ssnr(s(1:length(s_o3),1), s_o3, fs);
STOI3 = stoi(s(1:length(s_o3),1), s_o3, fs)
imp = SSNR3 - SSNR
%% Sound
sound(s_o3, fs)