clear all;
close all;

s = audioread('AudioFiles/clean_speech.wav');
n = 0.5*audioread('AudioFiles/babble_noise.wav');
ind = 1:70000;
fs = 16000;
l = 15;
o = 60;
win = 4;
P_l = 0;
% snr_ip = snr(s(ind), n(ind))

y = s(ind)+ n(ind);
% snr_ip1 = snr(y, fs)

Y = stft(y, win, l, o, 1, fs);
N = stft(n(ind), win, l, o, 1, fs);

Pyy = Bartlett_P(Y, 20);
Pnn = Bartlett_P(N, 20);

fixed_apriori_snr_db =15;
fixed_apriori_snr = 10^(fixed_apriori_snr_db/10);

%Initialize the noise estimate by assuming first eight frames are noise
varw_hat = zeros(size(Pyy));
varw_hat(:,1) =  Pyy(:,1);

for i = 2:size(Pyy,2)
p_h1_giveny = 1./(1+(1+fixed_apriori_snr)*exp((-Pyy(:,i)*fixed_apriori_snr)./(varw_hat(:,i-1)*(1+fixed_apriori_snr))));
p_alpha = 0.9;
P_l = p_alpha*P_l + (1-p_alpha)*p_h1_giveny;
if (P_l>0.99)
    p_h1_giveny = min(p_h1_giveny,0.99);
end
mmse_n= (1-p_h1_giveny).*Pyy(:,i)+p_h1_giveny.*varw_hat(:,i-1);
alpha = 0.8;
varw_hat(:,i) = alpha*varw_hat(:,i-1)+(1-alpha)*mmse_n;
end

true_varw = zeros(size(Pnn,2),1);
true_varw(1) = mean(Pnn(:,1));
for i = 2:size(Pnn,2)
    true_varw(i) = 0.9*mean(Pnn(:,i-1))+0.1*mean(Pnn(:,i));
end
figure
t = linspace(0,10800,size(Pnn,2))/100;
plot(t,10*log10(true_varw))
hold on
plot(t,10*log10(mean(varw_hat)))
xlabel('time (s)')
ylabel('\sigma^2_w (db)')
legend('true noise variance','SPP estimated noise variance')

S_w = Wiener(Pyy, varw_hat, Y);
s_out2 = stift(S_w, win, l, o, 1, fs);
sound(s_out2, fs)
% S_ps = Spectral_Subtraction(Pyy, varw_hat, Y);
% s_out1 = stift(S_ps, win, l, o, 1, fs);
% sound(s_out1, fs)
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