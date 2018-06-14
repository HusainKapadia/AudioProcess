function n = genNoise(type, dB, len)
    switch type
        case 1
            n = wgn(len, 1, dB);
        case 2
            n = 0.5*audioread('AudioFiles/babble_noise.wav');
        case 3
            n = 50*audioread('AudioFiles/aritificial_nonstat_noise.wav');
        case 4
            n = audioread('AudioFiles/Speech_shaped_noise.wav');
    end
end