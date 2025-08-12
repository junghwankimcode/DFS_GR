function pslr_dB = computePSLR(seq)

    % Assume frequency-domain input → convert to time-domain
    seq_time = (seq);

    % Autocorrelation
    [acf, lags] = xcorr(seq_time);
    acf= acf./max(abs(acf));
    % Main peak (zero lag)
    centerIdx = find(lags == 0, 1);
    main_peak = abs(acf(centerIdx));

    % Remove the main peak for sidelobe calculation
    sidelobes = [acf(1:centerIdx-1); acf(centerIdx+1:end)];

    % Max sidelobe magnitude (scalar)
    max_sidelobe = max(abs(sidelobes(:)));

    % Compute PSLR (in dB) ← note: now it's max_sidelobe / main_peak
    if main_peak == 0
        pslr_dB = Inf; % handle edge case
    else
        pslr_dB = 20 * log10(max_sidelobe / main_peak);
    end
end
