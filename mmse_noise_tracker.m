function [noise_psd, clean_fft] = mmse_noise_tracker(Pyy, Y)
noise_psd(:,1) = mean(Pyy(:,1:5),2);
smoothed_varw_hat = zeros(size(Pyy));
time_frame = 103; %safety frame 0.8 seconds
alpha = 0.98;
beta = 0.8;
S_hat = 0;

for i = 1:size(Pyy,2)
    if(i==1)
        apost_snr = (Pyy(:,i))./noise_psd(:,i);
    else
        apost_snr = (Pyy(:,i))./noise_psd(:,i-1);
    end
    apriori_ml = max(apost_snr-1,eps);
    mmse_n = (1./(1+apriori_ml).^2 + apriori_ml./((1+apriori_ml).*apost_snr)).*Pyy(:,i);
    if(i ~= 1)
        apriori_dd = max(alpha*abs(S_hat).^2./noise_psd(:,i-1) + (1-alpha)*apriori_ml,eps);
        G = gammainc(1./(1+apriori_dd),2);
        B = (1+apriori_dd).*G + exp(-1./(1+apriori_dd));
        noise_psd(:,i) = mmse_n./B;
        smoothed_varw_hat(:,i) = beta*smoothed_varw_hat(:,i-1)+(1-beta)*noise_psd(:,i);
    end
    
    %Safety Net
    if i>time_frame
        noise_psd(:,i) = max(smoothed_varw_hat(:,i),min(Pyy(:,i-time_frame:i),[],2));
    else
        noise_psd(:,i) = max(smoothed_varw_hat(:,i),min(Pyy(:,1:i),[],2));
    end
    apost_snr = (Pyy(:,i))./noise_psd(:,i);
    apriori_ml = max(apost_snr-1,eps);
    if(i==1)
        apriori_dd = apriori_ml;
    else
        apriori_dd = max(alpha*abs(S_hat).^2./noise_psd(:,i-1) + (1-alpha)*apriori_ml,eps);
    end
    
    S_hat = Wiener2(apriori_dd, Y(:,i));
    clean_fft(:,i) = S_hat;
end
    
end