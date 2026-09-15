clear; clc;


nDrop = 1e2;          
nInit = 10;            % I: NumInitialization
maxIter = 10;          % T: NumIteration
K = 3;                % K: NumPaths
P1 = 64; P2 = 64;     % Grid Size
M = 4;                % 4-QAM (2 bits/symbol)
bitsPerSymbol = log2(M);
bitsPerGrid = P1 * P2 * bitsPerSymbol; 

% Proposed
phi = 1/13; 
C = exp(1i*2*pi*phi);
dSet = C.^[0 1 3 9]; 

EbN0vec = 20;
BER_proposed = zeros(length(EbN0vec), 1);
BER_perfectCSIR = zeros(length(EbN0vec), 1);
BER_10CSIR = zeros(length(EbN0vec),1);
BER_15CSIR = zeros(length(EbN0vec),1);
BER_20CSIR = zeros(length(EbN0vec),1);
BER_OTFS = zeros(length(EbN0vec),1);


toBits = @(syms) int2bit(syms(:),bitsPerSymbol);

for indexSNR = 1:length(EbN0vec)
    SNR = EbN0vec(indexSNR) + 10*log10(bitsPerSymbol);
    fprintf('Simulating Eb/N0 = %d dB...\n', EbN0vec(indexSNR));
    
    matchingCount = zeros(nDrop, 1);
    err_perfect = zeros(nDrop, 1);
    err_10 = zeros(nDrop,1);
    err_15 = zeros(nDrop,1);
    err_20 = zeros(nDrop,1);
    err_otfs = zeros(nDrop,1);
    
    parfor indexDrop = 1:nDrop
        
        [~, ~, ~, ~, H_true] = gsGen_2D(EbN0vec(indexSNR), K, P1, P2);
        
        % Proposed & Baseline 
        symbols_prop = randi(M, [P1, P2]); % 1 to M
        X_prop = dSet(symbols_prop);
        
        symbols_qam = randi([0 M-1], [P1, P2]);
        X_qam = qammod(symbols_qam, M, 'UnitAveragePower', true);
        
        Y_prop = awgn(H_true .* X_prop, SNR, 'measured'); 
        
        
        Y_qam = awgn(H_true .* X_qam, SNR, 'measured'); 

        %% ------
        best_match_drop = 0;
        bits_true_prop = toBits(symbols_prop - 1); % 1x392
        %temp_symbols = zeros(P1, P2);

        for init = 1:nInit % I
            % 1. 초기 추정치 (Random sequence)
            curr_X_hat = dSet(randi(M, [P1, P2]));
            
            for t = 1:maxIter % T
                % 2. channel refining
                H_est = Y_prop ./ curr_X_hat;
                H_mid = zeros(P1, P2);
                for r = 1:P1
                    [~, ~, ~, H_hankel] = makeHankel(H_est(r, :));
                    [~, ~, ~, H_proj] = projectToLowRankHankel(H_hankel, K);
                    H_mid(r, :) = unwrapHankel(H_proj, P2); 
                end
                H_final = zeros(P1, P2);
                for c = 1:P2
                    [~, ~, ~, H_hankel] = makeHankel(H_mid(:, c).');
                    [~, ~, ~, H_proj] = projectToLowRankHankel(H_hankel, K);
                    H_final(:, c) = unwrapHankel(H_proj, P1).';
                end
                
               
                X_raw = Y_prop ./ H_final; 
                X_refined = X_raw;
             
                % 3. Quantize
                temp_symbols = zeros(P1, P2);
                for r = 1:P1
                    for c = 1:P2
                        [~, idx] = min(abs(X_refined(r, c) - dSet));
                        temp_symbols(r, c) = idx;
                        curr_X_hat(r, c) = dSet(idx);
                    end
                end
            end
            
            % Bit matching
            bits_est = toBits(temp_symbols - 1);
            match = sum(bits_est == bits_true_prop);
            if match > best_match_drop, best_match_drop = match; 
            end
        end
        matchingCount(indexDrop) = best_match_drop;

        %% ---Baseline ---
        demod = @(x) qamdemod(x, M, 'UnitAveragePower', true);
        bits_qam_true = toBits(symbols_qam); % 1x392
        
        % Perfect channel
        sym_perf = demod(Y_qam ./ H_true);
        bits_perf = toBits(sym_perf); 

        %SNR_linear = 10^(SNR/10);
        % 10% error
        H_10= awgn(H_true, 10, "measured");


        %X_10 = (conj(H_10).*Y_qam)./(abs(H_10).^2+(1/SNR_linear));
        sym_10 = demod(Y_qam./H_10);
        bits_10 = toBits(sym_10);

        % % 15% error
        % H_15 = awgn(H_true,8,"measured");
        % sym_15 = demod(Y_qam./H_15);
        % bits_15 = toBits(sym_15);
        % 
        % % 20% error
        % H_20 = awgn(H_true, 7.5,"measured");
        % sym_20 = demod(Y_qam./H_20);
        % bits_20 = toBits(sym_20);

        % --- [D] Baseline: OTFS ---
        K_otfs = P2; % Frequency 축 (Delay Bin)
        L_otfs = P1; % Time 축 (Doppler Bin)
        
        % 1. DD domain symbols generation (K x L)
        symbols_otfs = randi([0 M-1], [K_otfs, L_otfs]);
        X_otfs_dd = qammod(symbols_otfs, M, 'UnitAveragePower', true);
        
        % 2. OTFS modulation (DD -> TF)
        % S_vec -> (K*L) x 1 vector. reshape -> (K x L) Grid size (Rows: Freq (delay), Cols: Time)
        S_otfs_vec = OTFSmod(X_otfs_dd); 
        S_otfs_tf_grid = reshape(S_otfs_vec, [K_otfs, L_otfs]); 
        
        % 3. Channel pass
        % H_true -> (P1 x P2) (Time x Freq), So, S_otfs_tf_grid -> transpose(.')
        Y_otfs_tf = awgn((H_true .* S_otfs_tf_grid.'), SNR, 'measured');
        
        % 4. Equalization
        %X_hat_tf_grid = Y_otfs_tf ./ H_10;
        SNR_linear = 10^(SNR/10);
        X_hat_tf_grid = (conj(H_10).*Y_otfs_tf)./(abs(H_10).^2+(1/SNR_linear));
        
        % 5. OTFS demodulation (TF -> DD)
        % OTFSdemod는 (K*L) x 1 
        X_hat_tf_vec = X_hat_tf_grid.'; % (Freq x Time)
        X_hat_dd = OTFSdemod(X_hat_tf_vec(:), [K_otfs, L_otfs]);
        
        % 6. Bit calculate
        sym_otfs_dec = demod(X_hat_dd);
        
        bits_otfs_true = toBits(symbols_otfs);
        bits_otfs_est = toBits(sym_otfs_dec);
        err_otfs(indexDrop) = sum(bits_otfs_est ~= bits_otfs_true);
        

        err_perfect(indexDrop) = sum(bits_perf ~= bits_qam_true);
        err_10(indexDrop) = sum(bits_10 ~= bits_qam_true);
        % err_15(indexDrop) = sum(bits_15 ~= bits_qam_true);
        % err_20(indexDrop) = sum(bits_20 ~= bits_qam_true);
    end
    
    BER_proposed(indexSNR) = 1 - mean(matchingCount) / bitsPerGrid;
    BER_perfectCSIR(indexSNR) = sum(err_perfect) / (nDrop * bitsPerGrid);
    BER_10CSIR(indexSNR) = sum(err_10) / (nDrop * bitsPerGrid);
    % BER_15CSIR(indexSNR) = sum(err_15) / (nDrop * bitsPerGrid);
    % BER_20CSIR(indexSNR) = sum(err_20) / (nDrop * bitsPerGrid);
    BER_OTFS(indexSNR) = sum(err_otfs) / (nDrop * bitsPerGrid);
end

%%
figure(3)
semilogy(EbN0vec, BER_proposed, 'b-s', 'LineWidth', 1.5); hold on;
    semilogy(EbN0vec, BER_perfectCSIR, 'r-o','LineWidth',1.5);
    semilogy(EbN0vec, BER_10CSIR, 'c-h','LineWidth',1.5);
    % semilogy(EbN0vec, BER_15CSIR,'m-+','LineWidth',1.5);
    % semilogy(EbN0vec, BER_20CSIR, 'k-d','LineWidth',1.5);
    semilogy(EbN0vec, BER_OTFS, 'g-^','LineWidth',1.5);
    xticks(EbN0vec)
    grid on;
    xlabel('E_b/N_0 [dB]');
    ylabel('BER');
    % legend("Proposed method","OFDM (Perfect CSIR)",'OFDM (10% CSIR Error)', ...
    %        'OFDM (15% CSIR Error)','OFDM (20% CSIR Error)', ...
    %        'OTFS (10% CSIR Error)','Location','best');
    legend("Proposed method","OFDM (Perfect CSIR)",'OFDM (10% CSIR Error)', ...
           'OTFS (10% CSIR Error)','Location','best');
