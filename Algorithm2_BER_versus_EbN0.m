clear; clc;

% ==== Parameters ====
nDrop = 1e3;
nInit = 10;
maxIter = 10;
Kvec = 3;
Lvec = [512]; 
M = 4;
R = 1;
C = exp(1i*2*pi/13);
dSet = C.^[0 1 3 9];  
EbN0vec = 0:5:20;
delta_f = 15e3;
tau_max = 1e-6;

% ==== Main Loop ====
for indexL = 1:length(Lvec)
    L = Lvec(indexL);
    
    for indexK = 1:length(Kvec)
        K = Kvec(indexK);
        
        BER_alt1 = zeros(length(Kvec), length(EbN0vec));
        BER_noEstimation = zeros(length(EbN0vec), 1);
        BER_withEstimation1 = zeros(length(EbN0vec), 1);
        BER_withEstimation3 = zeros(length(EbN0vec), 1);
        BER_withEstimation4 = zeros(length(EbN0vec), 1);
        
        for indexSNR = 1:length(EbN0vec)
            SNR = EbN0vec(indexSNR) + 10*log10(R*log2(M));

            
            matchingCount = zeros(nDrop, 1);
            err_noEst = zeros(nDrop,1);
            err_est1 = zeros(nDrop,1);
            err_est3 = zeros(nDrop,1);
            err_est4 = zeros(nDrop,1);
            
            parfor indexDrop = 1:nDrop
                % === Ground Truth ===
                symbols = randi(M, [1, L]);
                symbols2 = randi([0 M-1], 1, L);
                d = dSet(symbols);
                modulatedSignal = qammod(symbols2, M, 'UnitAveragePower', true);

                delay_p = unifrnd(0, tau_max, [1, K]);
                IT = randn(1, K) + 1i*randn(1, K);
                CR = exp(-2i * pi * delay_p * delta_f);
                h = IT * transpose(CR).^([0:L-1]);
                y = awgn(h .* d, SNR, 'measured');

                % --- Alt V1: Hard Quantization ---
                [~,~,~,Y] = makeHankel(y);
                array2 = convertToBinaryVector(symbols - 1, M);

                best_match = -inf;
                for altInit = 1:nInit
                    d_alt = dSet(randi(M, [1, L]));
                    [~,~,~,D_alt] = makeHankel(d_alt);
                    for iter = 1:maxIter
                        A_hat = Y ./ D_alt;
                        [~,~,~,A_proj] = projectToLowRankHankel(A_hat, K);
                        D_est = Y ./ A_proj;
                        [m, n] = size(D_est);
                        L_reconstructed = m + n - 1;
                        d_est_vec = zeros(1, L_reconstructed);
                        for k = 1:L_reconstructed
                            [i, j] = find(bsxfun(@plus, (1:m)', (1:n)) == k+1);
                            d_est_vec(k) = mean(D_est(sub2ind([m, n], i, j)));
                        end
                        symbols_alt = zeros(1, L);
                        for idx = 1:L
                            [~, ind] = min(abs(d_est_vec(idx) - dSet));
                            symbols_alt(idx) = ind;
                        end
                        d_alt = dSet(symbols_alt);
                        [~,~,~,D_alt] = makeHankel(d_alt);
                    end
                    array1 = convertToBinaryVector(symbols_alt - 1, M);
                    match = sum(array1 == array2);
                    if match > best_match
                        best_match = match;
                    end
                end
                matchingCount(indexDrop) = best_match;

                % --- Baseline BER  ---
                received = awgn(h.*modulatedSignal,SNR,'measured');  
                true_h = h;
                demod = @(x) qamdemod(x, M, 'UnitAveragePower', true);

                % 1) Perfect channel
                equalized = received ./ true_h;
                err_noEst(indexDrop) = sum(symbols2 ~= demod(equalized));

                % 2) 10% error
                h1 = awgn(true_h, 10, 'measured');
                equalized1 = received ./ h1;
                err_est1(indexDrop) = sum(symbols2 ~= demod(equalized1));

                % 3) 15% error
                h3 = awgn(true_h, 8, 'measured');
                equalized3 = received ./ h3;
                err_est3(indexDrop) = sum(symbols2 ~= demod(equalized3));

                % 4) 20% error
                h4 = awgn(true_h, 7.5, 'measured');
                equalized4 = received ./ h4;
                err_est4(indexDrop) = sum(symbols2 ~= demod(equalized4));
            end

            % ===  BER  ===
            totalBits = nDrop * L * log2(M);
            BER_alt1(indexK, indexSNR) = 1 - mean(matchingCount) / (L * log2(M));
            BER_noEstimation(indexSNR) = sum(err_noEst) / totalBits;
            BER_withEstimation1(indexSNR) = sum(err_est1) / totalBits;
            BER_withEstimation3(indexSNR) = sum(err_est3) / totalBits;
            BER_withEstimation4(indexSNR) = sum(err_est4) / totalBits;
        end
    end
    
    %% ==== Plot ====
    figure(1); clf;
    semilogy(EbN0vec, BER_alt1, 'b-s', 'LineWidth', 1.5); hold on;
    semilogy(EbN0vec, BER_noEstimation, 'r-o','LineWidth',1.5);
    semilogy(EbN0vec, BER_withEstimation1, 'c-h','LineWidth',1.5);
    semilogy(EbN0vec, BER_withEstimation3,'m-+','LineWidth',1.5);
    semilogy(EbN0vec, BER_withEstimation4, 'k-d','LineWidth',1.5);
    xticks(EbN0vec)
    grid on;
    xlabel('E_b/N_0 [dB]');
    ylabel('BER');
    legend("Proposed method","Perfect CSIR",'Imperfect CSIR (10% error)', ...
           'Imperfect CSIR (15% error)','Imperfect CSIR (20% error)', ...
           'Location','best');
    %title(sprintf('BER Comparison (L = %d, K = %d)', L, K));
end
