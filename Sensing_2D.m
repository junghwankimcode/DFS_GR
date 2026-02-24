clear; %clc;


% 1) Simulation parameters
P1 = 512; P2 = 64;                 % P1: slow-time (Doppler bins), P2: freq (Delay bins)
SNR_dB = -6;
numIterations = 1e2;

c_speed      = physconst("lightspeed");
carrier_freq = 28e9;
delta_f      = 480e3;
T            = 1/delta_f;           % OFDM symbol duration

% Bin resolutions (radar round-trip assumed)
range_per_bin = c_speed / (2 * P2 * delta_f);
vel_per_bin   = c_speed / (2 * carrier_freq * P1 * T);

K = 8                               % number of paths for TF-channel generator

% Proposed symbol set (Golomb ruler)
phi  = 1/13;
C    = exp(1i*2*pi*phi);
dSet = C.^[0 1 3 9];


% Result: [R_true, V_true, Rhat1,Vhat1, Rhat2,Vhat2, Rhat3,Vhat3, Rhat4,Vhat4]
res = zeros(numIterations, 6);

fprintf('Testing SNR = %d dB...\n', SNR_dB);

%[nn, mm] = ndgrid(0:P1-1, 0:P2-1);

% 2) Monte-Carlo loop

for iter = 1:numIterations

    % --- TF-domain channel (your original generator) ---
    % H_true is a 2D exponential sum channel on (n,m)-grid.
    [~, ~, ~, ~, H_true, R_true, V_true] = gsGen_2D_sensing(SNR_dB, K, P1, P2);

    % (1) Proposed
    X_prop = dSet(randi([1 4], P1, P2));
    Y_p    = awgn(H_true .* X_prop, SNR_dB, 'measured');

    % % (2) ZC (across frequency axis)
    % numbers = primes(P2);
    % r_cand  = numbers(gcd(numbers, P2) == 1);
    % r       = r_cand(randi(length(r_cand)));
    % zc_base = zadoffChuSeq(r, P2);
    % X_zc    = repmat(zc_base.', P1, 1);
    % Y_zc    = awgn(H_true .* X_zc, SNR_dB, 'measured');
    % 
    % % (3) M-seq (across frequency axis)
    % % NOTE: depending on your MATLAB version, mlseq input may be "order" not "length".
    % % If your mlseq expects order, replace: m_base = mlseq(9);  then ensure length==511.
    % m_base = mlseq(P2);
    % X_m    = repmat(m_base.', P1, 1);
    % Y_m    = awgn(H_true .* X_m, SNR_dB, 'measured');

    % (4) OTFS (4-QAM) - "OTFS-consistent" DD-domain single-tap model + DD matched filter
    % Use the SAME R_true, V_true so the physical truth is consistent.
    K_otfs  = K;    
    K_dB    = 6;

    [Rhat_otfs, Vhat_otfs] = otfs_rd_peak_multipath( ...
    R_true, V_true, K_otfs, K_dB, ...
    P1, P2, delta_f, T, carrier_freq, c_speed, ...
    SNR_dB);

    % 3) RD-map peak estimator for (1)-(3) in TF model
   
    % Common processing: diff = Y .* conj(X)  -> IFFT over freq (range) -> FFT over slow-time (Doppler)
    [Rhat_prop, Vhat_prop] = rd_peak_from_tf(Y_p,  X_prop, P1, P2, range_per_bin, vel_per_bin);
    % [Rhat_zc,   Vhat_zc]   = rd_peak_from_tf(Y_zc, X_zc,   P1, P2, range_per_bin, vel_per_bin);
    % [Rhat_m,    Vhat_m]    = rd_peak_from_tf(Y_m,  X_m,    P1, P2, range_per_bin, vel_per_bin);

    % store
    res(iter, :) = [R_true, V_true, ...
                    Rhat_prop, Vhat_prop, ...
                    Rhat_otfs, Vhat_otfs];
end


%% 4) Print summary
methods = {'Proposed', 'OTFS (4-QAM)'};
all_R_real = res(:,1); all_V_real = res(:,2);

all_R_est = cell(1,2);
all_V_est = cell(1,2);
for m = 1:2
    all_R_est{m} = res(:, 2*m+1);
    all_V_est{m} = res(:, 2*m+2);
end

line_sep = repmat('=', [1, 85]);
dash_sep = repmat('-', [1, 85]);

fprintf('\n%s\n', line_sep);
fprintf('%-18s | %-30s | %-30s\n', 'Method', 'Range Error (m)', 'Velocity Error (m/s)');
fprintf('%s\n', dash_sep);

for m = 1:2
    r_err = abs(all_R_real - all_R_est{m});
    v_err = abs(all_V_real - all_V_est{m});

    m_r = mean(r_err); s_r = std(r_err);
    m_v = mean(v_err); s_v = std(v_err);

    fprintf('%-18s | Mean: %10.5f (Std: %10.5f) | Mean: %10.5f (Std: %10.5f)\n', ...
            methods{m}, m_r, s_r, m_v, s_v);
end
fprintf('%s\n', line_sep);


%% ============================================================
%                     Functions


function [Rhat, Vhat] = rd_peak_from_tf(Y, X, P1, P2, range_per_bin, vel_per_bin)
    % RD-map from TF-like grid:
    % diff = Y .* conj(X) ~ H + noise (if |X|=1)
    % range: IFFT along dim-2, doppler: FFT along dim-1
    diff = Y .* conj(X);
    RD   = fft(ifft(diff, P2, 2), P1, 1);

    % Peak search
    [~, idx] = max(abs(RD(:)).^2);
    [v_idx, r_idx] = ind2sub(size(RD), idx);

    % Convert indices to lags
    lag_r = mod(r_idx - 1, P2);

    half_P1 = floor(P1/2);
    lag_v = mod(v_idx - 1 + half_P1, P1) - half_P1;

    % Physical mapping
    Rhat = lag_r * range_per_bin;
    Vhat = lag_v * vel_per_bin;
end

%%
function [Rhat, Vhat, RD, info] = otfs_rd_peak_multipath( R_true, V_true, K, K_dB, ...
    P1, P2, delta_f, T, fc, c, ...
    SNR_dB)

% OTFS RD-map baseline (DD-domain):
% Y_DD = sum_p alpha_p * circshift(X_DD, [k_p, ell_p]) + W
% RD   = ifft2( fft2(Y_DD) .* conj(fft2(X_DD)) )
% Peak -> (k,ell) -> (v,R)

    % 0) constants
    range_per_bin = c / (2 * P2 * delta_f);
    vel_per_bin   = c / (2 * fc * P1 * T);
    halfP1 = floor(P1/2);

    % 1) DD symbols (4-QAM)
    sym  = randi([0 3], [P1, P2]);
    X_DD = qammod(sym, 4, 'UnitAveragePower', true);

    % 2) multipath parameters
    % physical ranges/velocities: LoS is (R_true,V_true), others random
    R_limit = 200;
    V_limit = 20;
    R_all = [R_true, unifrnd(0, R_limit, [1, max(K-1,0)])];
    V_all = [V_true, unifrnd(-V_limit, V_limit, [1, max(K-1,0)])];

    % Rician power split
    K_lin = 10^(K_dB/10);

    % allocate
    Y_DD_ideal = zeros(P1, P2);
    k_bins = zeros(1,K);
    l_bins = zeros(1,K);
    alpha  = zeros(1,K);

    for p = 1:K
        % --- gains (Rician) ---
        if p == 1
            % deterministic dominant (unit-magnitude with random phase) scaled by sqrt(K/(K+1))
            alpha(p) = sqrt(K_lin/(K_lin+1)) * exp(1j*2*pi*rand);
        else
            % Rayleigh scatter part total power = 1/(K+1), equally split
            alpha(p) = sqrt(1/(K_lin+1)) * (randn + 1j*randn) / sqrt(2*(K-1));
        end

        % --- physical -> bins (GRID-ALIGNED) ---
        tau = 2 * R_all(p) / c;          % round-trip delay
        fD  = 2 * V_all(p) * fc / c;     % round-trip Doppler

        ell = mod(round(tau * delta_f * P2), P2);   % delay bin 0..P2-1

        k_raw = round(fD * T * P1);                 % possibly negative
        k_cent = mod(k_raw + halfP1, P1) - halfP1;  % centered bin in [-P1/2, P1/2)
        k_shift = mod(k_cent, P1);                  % shift index 0..P1-1

        k_bins(p) = k_cent;     % store centered index for interpretation
        l_bins(p) = ell;

        % --- DD convolution: circshift ---
        Y_DD_ideal = Y_DD_ideal + alpha(p) * circshift(X_DD, [k_shift, ell]);
    end

    % 3) noise
    Y_DD = awgn(Y_DD_ideal, SNR_dB, 'measured');

    % 4) RD-map (DD matched filter)
    RD = ifft2( fft2(Y_DD) .* conj(fft2(X_DD)) );

    % 5) strongest peak detection
    [~, idx] = max(abs(RD(:)).^2);
    [k_idx, l_idx] = ind2sub(size(RD), idx);

    % lags
    lag_l = mod(l_idx - 1, P2);
    lag_k0 = mod(k_idx - 1, P1);
    lag_k  = mod(lag_k0 + halfP1, P1) - halfP1;

    % 6) bin -> physical
    Rhat = lag_l * range_per_bin;
    Vhat = lag_k * vel_per_bin;

    % 7) debug info
    info.range_per_bin = range_per_bin;
    info.vel_per_bin   = vel_per_bin;
    info.true_R_all = R_all;
    info.true_V_all = V_all;
    info.true_k_bins = k_bins;
    info.true_l_bins = l_bins;
    info.true_alpha  = alpha;
    info.peak_k = lag_k;
    info.peak_l = lag_l;
end
