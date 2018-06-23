function [noise_psd] = mmse_noise_tracker_spp(Pyy)
%Initialize the noise estimate by assuming first eight frames are noise
noise_psd = zeros(size(Pyy));
noise_psd(:,1) =  mean(Pyy(:,1:5),2);
P_l = 0;

fixed_apriori_snr_db =15;
fixed_apriori_snr = 10^(fixed_apriori_snr_db/10);

for i = 2:size(Pyy,2)
p_h1_giveny = 1./(1+(1+fixed_apriori_snr)*exp((-Pyy(:,i)*fixed_apriori_snr)./(noise_psd(:,i-1)*(1+fixed_apriori_snr))));
p_alpha = 0.9;
P_l = p_alpha*P_l + (1-p_alpha)*p_h1_giveny;
if (P_l>0.99)
    p_h1_giveny = min(p_h1_giveny,0.99);
end
mmse_n= (1-p_h1_giveny).*Pyy(:,i)+p_h1_giveny.*noise_psd(:,i-1);
alpha = 0.85;
noise_psd(:,i) = alpha*noise_psd(:,i-1)+(1-alpha)*mmse_n;
end
end
