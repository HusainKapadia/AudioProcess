function lrt = vad(apriori_snr,apost_snr)

    lambda_k = (1./(1+apriori_snr)).*exp((apost_snr.*apriori_snr./(1+apriori_snr)));
    lrt = mean(log(lambda_k));

end