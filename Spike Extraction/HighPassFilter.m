function w = HighPassFilter(wave)

    load Filter.mat

    wave = filter(Hhp5,wave);

    l = length(wave);

    wave_fft = fft(wave);

    mean_wave = filter(ones(400,1),400,abs(wave_fft));

    k = 0;
    for i = 2:(l/2 - 5)
        h = abs(wave_fft(i));
        if h > 2.5 * mean_wave(i)
            wave_fft(i) = 0;
            wave_fft(l - i + 2) = 0;
        end
    end

    %figure; plot(abs(wave_fft));

    w = ifft(wave_fft,'symmetric');