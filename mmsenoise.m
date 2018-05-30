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
S = stft(s(ind), win, l, o, 1, fs);
Pyy = Bartlett_P(Y, 8);
Pnn = Bartlett_P(N, 8);
Pss = Bartlett_P(S,8);
%Initialize the noise estimate by assuming first eight frames are noise
varw_hat = zeros(size(Pyy));
varw_hat =  Pyy(:,1);
Pnn(:,1) = varw_hat;
ebs = zeros(size(Pyy));
smoothed_varw_hat = zeros(size(Pyy));
%Update noise using safety net
time_interval = 0.8; %seconds
time_frame = floor(0.8/0.015); %safety frame
Rprior=-40:1:100;% dB
[tabel_inc_gamma ]=tab_inc_gamma(Rprior,2);
for i = 2:720
    ebs = max((Pyy(:,i))./varw_hat(:,i-1)-1,0);
    aposteriori_SNR = Pyy(:,i)./varw_hat(:,i-1);
    mmse_n = (1./(1+ebs).^2+ebs./((1+ebs).*aposteriori_SNR)).*Pyy(:,i);
    alpha = 0.8;
    ebs_dd = max(alpha*varw_hat(:,i-1)/mean(varw_hat(:,i-1)) + (1-alpha)*abs(Pyy(:,i))/mean(varw_hat(:,i-1))-1,0);
    
    G = gammainc(2,1./(1+ebs_dd),'scaledupper');
    B = (1+ebs_dd).*G+exp(-1./(1+ebs_dd));
    varw_hat(:,i) = mmse_n.*B;
    beta = 0.8;
    smoothed_varw_hat(:,i) = beta*varw_hat(:,i-1)+(1-beta)*varw_hat(:,i);
    if i>time_frame
    varw_hat(:,i) = max(smoothed_varw_hat(:,i),min(Pyy(:,i-time_frame:i),[],2));
    end
    
end
true_varw = zeros(size(Pnn,2),1);
true_varw(1) = mean(Pnn(:,1));
for i = 2:size(Pnn,2)
    true_varw(i) = 0.9*mean(Pnn(:,i-1))+0.1*mean(Pnn(:,i));
end
figure
t = linspace(0,10800,720)/100;
plot(t,10*log10(true_varw))
hold on
plot(t,10*log10(mean(varw_hat)))
%%
aposteriori_SNR = (Pyy./Pnn);
ebs = zeros(size(Pyy));
for i = 2:720
    ebs(:,i) = max(abs(Pyy(:,i))./Pnn(:,i-1)-1,0);
end
mmse_n = (1./(1+ebs).^2+ebs./((1+ebs).*aposteriori_SNR)).*Pyy;
alpha = 0.8;
for i = 2:720
ebs_dd(:,i) = alpha*Pnn(:,i-1)/var(N(:,i)) + (1-alpha)*max((abs(Pyy(:,i))/var(N(:,i)))-1,0);
end
G = tab_inc_gamma(2,1./(1+ebs_dd));
B = (1+ebs_dd)*G+exp(-1./(1+ebs_dd));
varw_hat = mmse_n.*B;
beta = 0.8;
for i = 2:720
smoothed_varw_hat(:,i) = beta*varw_hat(:,i-1)+(1-beta)*varw_hat(:,i);
end
%Update noise using safety net
time_interval = 0.8; %seconds
time_frame = floor(0.8/0.015); %safety frame
for i = 1:size(Pyy,2) 
    varw_hat(:,i) = max(smoothed_varw_hat(:,i),min(Pyy(:,i:i+time_frame),[],2));
end