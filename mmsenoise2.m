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
% Pnn(:,1) = varw_hat;
%ebs = zeros(size(Pyy));
smoothed_varw_hat = zeros(size(Pyy));
%Update noise using safety net
time_interval = 0.8; %seconds
time_frame = 103; %safety frame
%varw_hat(:,floor(1600/l)) = mean(Pyy(:,1:floor(1600/l)),2);
alpha = 0.98;
beta = 0.8;
[~, apriori_snr] = SNR_estimates(Pyy(:,floor(1600/l)), varw_hat(:,floor(1600/l)-1), 'ML');
S_hat = 0;
for i = 1:size(Pyy,2)
    if(i==1)
        apost_snr = (Pyy(:,i))./varw_hat(:,i);
    else
        apost_snr = (Pyy(:,i))./varw_hat(:,i-1);
    end
    ebs = max(apost_snr-1,eps);
    mmse_n = (1./(1+ebs).^2 + ebs./((1+ebs).*apost_snr)).*Pyy(:,i);
    if(i == 1)
        ebs_dd = ebs;    
    else
        ebs_dd = max(alpha*abs(S_hat).^2./varw_hat(:,i-1) + (1-alpha)*ebs,eps);
    end
    G = gammainc(1./(1+ebs_dd),2);
    B = (1+ebs_dd).*G + exp(-1./(1+ebs_dd));
    varw_hat(:,i) = mmse_n./B;
    if(i ~= 1)
        smoothed_varw_hat(:,i) = beta*varw_hat(:,i-1)+(1-beta)*varw_hat(:,i);
    end

    %Safety Net
    if i>time_frame
        varw_hat(:,i) = max(smoothed_varw_hat(:,i),min(Pyy(:,i-time_frame:i),[],2));
    end
    %S_hat = ebs_dd.*varw_hat(:,i-1);
    apost_snr = (Pyy(:,i))./varw_hat(:,i);
    ebs = max(apost_snr-1,eps);
    if(i==1)
        ebs_dd = ebs;
    else
        ebs_dd = max(alpha*abs(S_hat).^2./varw_hat(:,i-1) + (1-alpha)*ebs,eps);
    end
    
    S_hat = Wiener2(ebs_dd, Y(:,i));
    S_w(:,i) = S_hat;
end
true_varw = zeros(size(Pnn,2),1);
true_varw(1) = mean(Pnn(:,1));
for i = 2:size(Pnn,2)
    true_varw(i) = beta*mean(Pnn(:,i-1))+(1-beta)*mean(Pnn(:,i));
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

%S_w = Wiener1(Pyy, varw_hat, Y);
s_out2 = stift(S_w, win, l, o, 1, fs);

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