function [noise_psd,clean_psd] = vad_noise_tracker(Pyy, Y)
%silent_frames = zeros(l,1);
%Initialize the noise estimate by assuming first eight frames are noise

noise_psd = zeros(size(Pyy));
noise_psd(:,1) = Pyy(:,1);
alpha = 0.99;
beta = 0.9;
S_hat = 0;
%lrt_expanded = zeros(l,1);
for i = 2:size(Pyy,2)
    
    if(i==1)
        apost_snr = (Pyy(:,i))./noise_psd(:,i);
    else
        apost_snr = (Pyy(:,i))./noise_psd(:,i-1);
    end
    apriori_ml = max(apost_snr-1,eps);
    if(i ~= 1)
        apriori_dd = max(alpha*abs(S_hat).^2./noise_psd(:,i-1) + (1-alpha)*apriori_ml,eps);
    else
        apriori_dd = apriori_ml;
    end
    lrt(i) = vad(apriori_dd,apost_snr);
    if(isinf(lrt(i)))
        lrt(i) = lrt(i-1);
    end  
    if (lrt(i)>0.04)
        noise_psd(:,i) = noise_psd(:,i-1);
        %silent_frames((i-1)* (1 - 0.01*o)*frame+1:frame+(i-1)* (1 - 0.01*o)*frame) = ones(frame,1);
    else
        noise_psd(:,i) = beta*noise_psd(:,i-1)+(1-beta)*Pyy(:,i);
    end
    S_hat = Wiener2(apriori_dd,Y(:,i));
    clean_psd(:,i) = S_hat;
    %lrt_expanded((i-1)* (1 - 0.01*o)*frame+1:frame+(i-1)* (1 - 0.01*o)*frame) = ones(frame,1)*lrt(i);
end
%silent_frames(silent_frames==1) = nan;

end
