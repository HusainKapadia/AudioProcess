function n = genNoise(type, dB, len)
    switch type
        case 1
            n = wgn(len, 1, -dB);
        case 2
            n = dB*audioread('AudioFiles/babble_noise.wav')./50;
        case 3
            n = dB*audioread('AudioFiles/aritificial_nonstat_noise.wav');
        case 4
            n = dB*audioread('AudioFiles/Speech_shaped_noise.wav')./100;
    end
end