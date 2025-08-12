function snr_dB = getSNRfromError(relativeErrorPercent)
% getSNRfromError  Compute required SNR [dB] for a given relative error (%)
%
%   Input:
%       relativeErrorPercent - Target relative error in percentage (e.g., 10 for 10%)
%   Output:
%       snr_dB - Required SNR in dB to achieve the given relative error

    epsilon = relativeErrorPercent / 100;  % Convert percent to ratio
    snr_linear = 1 / (epsilon^2);          % SNR = 1 / error^2
    snr_dB = 10 * log10(snr_linear);       % Convert to dB
end
