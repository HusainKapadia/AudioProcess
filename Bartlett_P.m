function P = Bartlett_P(Y, M)
    P = zeros(size(Y));
    m = mod(size(Y,2),M);

    for i = 0:M:size(Y,2)
        if(i>(size(Y,2)-M))
            if(m == 0)
                P(:,i) = sum(abs(Y(:,i)).^2, 2);
            else 
                P(:,i+1:i+m) = repmat(sum(abs(Y(:,i+1:i+m)).^2, 2)/m, 1, m);
            end
        else
            P(:,i+1:i+M) = repmat(sum(abs(Y(:,i+1:i+M)).^2, 2)/M, 1, M);
        end
    end
    
end