function Y = stft(x, win, frame_len, overlap, mic, fs)
     L = frame_len*fs/1000;                %frame length                         %percent overlap
     D = (1 - 0.01*overlap)*L;             %start index for overlap
     K = round(1 + floor((length(x)-L)/D));%number of sections
     Y = zeros(L, K, mic);
     
     switch win
         case 1 
             w = ones(L,1);
         case 2 
             w = hamming(L);
         case 3 
             w = hanning(L);
         case 4 
             w = bartlett(L);
         case 5 
             w = blackman(L);
     end

     for j = 1:mic
         n1 = 1;                           %start index
         for i=1:K
             xw(:, j) = x(n1:n1+L-1, j).*w/norm(w);
             Y(:, i, j) = fft(xw(:, j), L);
             n1 = n1 + D;
         end
     end

end