function pslr_dB = computePSLR_circular(seq)

   
    % Circular autocorrelation은 아래처럼 계산
    seq=(seq);
    acf_circ = ifft(abs(fft(seq)).^2);  % circular autocorrelation
    acf_circ=acf_circ./max(abs(acf_circ));
    main_peak = abs(acf_circ(1));       % lag = 0에 해당
    sidelobes = abs(acf_circ(2:end));   % 나머지는 sidelobe

    max_sidelobe = max(sidelobes);

    if main_peak == 0
        pslr_dB = Inf;
    else
        pslr_dB = 20 * log10(max_sidelobe / main_peak); % sidelobe / mainlobe
    end
end
