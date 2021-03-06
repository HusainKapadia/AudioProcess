clear all;
close all;

s = audioread('AudioFiles/clean_speech.wav');
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
<<<<<<< HEAD

varw_hat(:,1) =  mean(Pyy(:,1:(floor(1/0.015))),2);
=======
varw_hat =  Pyy(:,1:floor(1600/l));
% Pnn(:,1) = varw_hat;
%ebs = zeros(size(Pyy));
>>>>>>> 86e6620a1137b52fb45f93af9039ace8bd842c57
smoothed_varw_hat = zeros(size(Pyy));
%Update noise using safety net
time_interval = 0.8; %seconds
time_frame = floor(0.8/0.015);

<<<<<<< HEAD
for i = 2:size(Pyy,2)
    %S_w = Wiener(Pyy(:,i-1),varw_hat(:,i-1), Y(:,i));
    aposteriori_SNR = Pyy(:,i)./varw_hat(:,i-1);
    %ALPHA = 0.8;
    %apriori_snr=max(ALPHA*( abs(S_w).^2./varw_hat(:,i-1))+(1-ALPHA)*(aposteriori_SNR -1),0);
    
    apriori_snr = max((Pyy(:,i))./varw_hat(:,i-1)-1,0);
    
    mmse_n = (1./(1+apriori_snr).^2 + apriori_snr./((1+apriori_snr).*aposteriori_SNR)).*Pyy(:,i);
    alpha = 0.9;
=======
for i = floor(1600/l):size(Pyy,2)
    apripori_snr = max((Pyy(:,i))./varw_hat(:,i-1)-1,0);
    aposteriori_SNR = Pyy(:,i)./varw_hat(:,i-1);
    mmse_n = (1./(1+apripori_snr).^2 + apripori_snr./((1+apripori_snr).*aposteriori_SNR)).*Pyy(:,i);
    alpha = 0.8;
>>>>>>> 86e6620a1137b52fb45f93af9039ace8bd842c57
    ebs_dd = max(alpha*varw_hat(:,i-1)/mean(varw_hat(:,i-1)) + (1-alpha)*abs(Pyy(:,i))/mean(varw_hat(:,i-1))-1,0);
    
    G = gammainc(2,1./(1+ebs_dd),'scaledupper');
    B = (1+ebs_dd).*G+exp(-1./(1+ebs_dd));
    varw_hat(:,i) = mmse_n.*B;
    beta = 0.8;
    smoothed_varw_hat(:,i) = beta*varw_hat(:,i-1)+(1-beta)*varw_hat(:,i);
    %Safety Net
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
t = linspace(0,10800,727)/100;
plot(t,10*log10(true_varw))
hold on
plot(t,10*log10(mean(varw_hat)))
xlabel('time (s)')
ylabel('\sigma^2_w (db)')
legend('true noise variance','mmse estimated noise variance')
%%
S_ps = Spectral_Subtraction(Pyy, varw_hat, Y);
s_out1 = stift(S_ps, win, l, o, 1, fs);

<<<<<<< HEAD
S_w = Wiener(Pyy, varw_hat, Y);
s_out2 = stift(S_w, win, l, o, 1, fs);
=======
S_w = Wiener1(Pyy, varw_hat, Y);
s_out2 = stift(S_ps, win, l, o, 1, fs);
>>>>>>> 86e6620a1137b52fb45f93af9039ace8bd842c57

%sound(s_out1, fs)
sound(s_out2, fs)
figure
subplot(3,1,1)
plot(s_out2);
axis([0 70000 -0.5 0.5])
title('Estimated Clean Speech 1')
subplot(3,1,2)
plot(y);
axis([0 70000 -0.5 0.5])
title('Noisy Speech')
subplot(3,1,3)
plot(s);
axis([0 70000 -0.5 0.5])
title('Clean Speech')