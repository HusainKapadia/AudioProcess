function [apost_SNR, aprior_SNR] = SNR_estimates(Pyy, Pnn, init, S_prev, a)
    
    apost_SNR = Pyy./Pnn;
    
    if init=='ML' & nargin == 3
        aprior_SNR = max(apost_SNR-1,0);
    elseif init =='DD'    
        aprior_SNR = a*abs(S_prev).^2./Pnn + (1-a)*max(Pyy./Pnn-1, 0);
    end
    
end