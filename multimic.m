clc;
close all;

s = audioread('AudioFiles/clean_speech.wav');
type = 1;
db = 20;
n1 = genNoise(type, db, length(s));
n2 = genNoise(type, db, length(s));
n3 = genNoise(type, db, length(s));
n4 = genNoise(type, db, length(s));

fs = 16000;                     %sampling frequency
ind = 1:96240;                  %length of speech
m = 4;                          %number of mics
l = 15;                         %frame length in ms
o = 60;                         %percent overlap

%% Generate Noisy Speech
c = 340;                        %speed of sound (m/s)
a = 40;                         %target direction angle (degrees)
x = 0.1*(1:m-1);                       %distance between mics (m) 
tau = x.*fs.*cos(pi*a/180)./c;     %Time delay (s)

s = repmat(s(ind), 1, m);          %clean speech
n = [n1(ind) n2(ind) n3(ind) n4(ind)];             %noise
    
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
%[noise_psd] = mmse_noise_tracker(Pse, S_e);
%[~, ~, noise_psd] = vad_noise_tracker(Pse, S_e, size(s(:,1)));
[noise_psd, S_o] = mmse_noise_tracker_spp(Pse,Y,'Wiener');

%%STIFT with overlap add
s_o = stift(S_o, 4, l, o, 1, fs);

%% Plots
t = (1:96240)/fs;
figure()
subplot(3,1,1)
plot(t,y(ind)), title('Noisy Speech');
axis([1 96240/fs -0.5 0.5]);

subplot(3,1,2)
plot(t,s_o(ind)), title('Corrected Speech');
axis([1 96240/fs -0.5 0.5]);

subplot(3,1,3)
plot(t,s(ind)), title('Clean Speech');
axis([1 96240/fs -0.5 0.5]);
%% Sound
sound(s_o, fs)