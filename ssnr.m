function s = ssnr(clean,estimated,fs)
     L = 15*fs/1000;                       % 15 ms frame length
     M = floor(size(estimated, 1)/L);      %number of sections
     s = 0;
     for i = 1:M-L
        s = s + min(max(10*log(sum(clean(i*L:(i+1)*L).^2)./...
            sum(estimated(i*L:(i+1)*L).^2)),-10), 35);
     end
     s = s/M;
end
    