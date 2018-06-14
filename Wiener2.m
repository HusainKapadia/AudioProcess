function S = Wiener2(SNR, Y)
    S = zeros(size(Y));
    
    for i = 1:size(Y,2)
        S(:,i) = (SNR./(1+SNR)).*Y(:,i);
    end
    
end