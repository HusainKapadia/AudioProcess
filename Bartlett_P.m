function P = Bartlett_P(Y, M)
    P = zeros(size(Y));
    
    for i = 0:M:(size(Y,2)-M)
        P(:,i+1:i+M) = repmat(sum(abs(Y(:,i+1:i+M)).^2, 2)/M, 1, M);
    end
end