function S = Wiener1(Pyy, Pnn, Y)
    S = zeros(size(Y));
    
    for i = 1:size(Y,2)
        S(:,i) = (1-Pnn(:,i)./Pyy(:,i)).*Y(:,i);
    end
    
end