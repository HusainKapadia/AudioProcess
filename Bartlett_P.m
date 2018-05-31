function P = Bartlett_P(Y, M)
    P = zeros(size(Y));
    m = mod(size(Y,2),M);
    for i = 0:M:(size(Y,2)-M)
        P(:,i+1:i+M) = repmat(sum(abs(Y(:,i+1:i+M)).^2, 2)/M, 1, M);
        if(i==(size(Y,2)-M))
            P(:,i+1:i+M-m) = repmat(sum(abs(Y(:,i+1:i+M-m)).^2, 2)/(M-m), 1, M-m);
        end
    end
end