function S = Spectral_Subtraction(Pyy, Pnn, Y)
    
    S_mag = zeros(size(Y));
    phase = angle(Y);
    
    for i = 1:size(Y,2)
        S_mag(:,i) = sqrt(max(1-(Pnn(:,i)./Pyy(:,i)), 0)).*abs(Y(:,i));
    end
    
    S = S_mag.*exp(1i*phase);
end