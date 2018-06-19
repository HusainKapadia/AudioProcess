function [noise_psd,clean_psd] = mmse(Y)
M = 8;
Pyy = Bartlett_P(Y, M);
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
    ebs = max(apost_snr-1,eps);
    mmse_n = (1./(1+ebs).^2 + ebs./((1+ebs).*apost_snr)).*Pyy(:,i);
    if(i == 1)
        ebs_dd = ebs;
    else
        ebs_dd = max(alpha*abs(S_hat).^2./noise_psd(:,i-1) + (1-alpha)*ebs,eps);
        G = gammainc(1./(1+ebs_dd),2);
        B = (1+ebs_dd).*G + exp(-1./(1+ebs_dd));
        noise_psd(:,i) = mmse_n./B;
    end
    
    if(i ~= 1)
        smoothed_varw_hat(:,i) = beta*smoothed_varw_hat(:,i-1)+(1-beta)*noise_psd(:,i);
    end
    
    %Safety Net
    if i>time_frame
        noise_psd(:,i) = max(smoothed_varw_hat(:,i),min(Pyy(:,i-time_frame:i),[],2));
    end
    %S_hat = ebs_dd.*varw_hat(:,i-1);
    apost_snr = (Pyy(:,i))./noise_psd(:,i);
    ebs = max(apost_snr-1,eps);
    if(i==1)
        ebs_dd = ebs;
    else
        ebs_dd = max(alpha*abs(S_hat).^2./noise_psd(:,i-1) + (1-alpha)*ebs,eps);
    end
    
    S_hat = Wiener2(ebs_dd, Y(:,i));
    clean_psd(:,i) = S_hat;
end
