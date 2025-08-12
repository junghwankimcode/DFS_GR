clear; clc;
%% Setting the Proposed Method (Hankel Decoder)
m = 0.5:0.01:2;
m = m(67);
phi = 0.01:0.01:1;
phi = 1/13;%phi(41);
%C = m*exp(1i*2*pi*(phi));
C = exp(1i*2*pi*(phi));
%dSet = C.^[0 1 4 6];
dSet = C.^[0 1 3 9];
delta_f = 15*1e3;  
tau_max = 1e-6;
%tau_max = 1/delta_f;

%% M-QAM BER Simulator with Channel Estimation
% Configuration Parameters
M = 4; % Modulation order 
EbN0_dB = 0:5:25; % Range of Eb/N0 values in dB
numIterations = 1e3; % Number of iterations

L = 8; % Length of received signal per iteration
K = 2;

BER_noEstimation = zeros(length(EbN0_dB), 1); % Without channel estimation
BER_withEstimation1 = zeros(length(EbN0_dB), 1); % With channel estimation
BER_withEstimation2 = zeros(length(EbN0_dB), 1); % With channel estimation
BER_withEstimation3 = zeros(length(EbN0_dB), 1); % With channel estimation
BER_withEstimation4 = zeros(length(EbN0_dB), 1); % With channel estimation

% Simulation Loop over Eb/N0
for idx = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(idx) / 10); % Linear Eb/N0
    noisePower = 1 / (2 * log2(M) * EbN0); % Noise power per dimension
    numErrors_noEstimation = 0;
    numErrors_withEstimation1 = 0;
    numErrors_withEstimation2 = 0;
    numErrors_withEstimation3 = 0;
    numErrors_withEstimation4 = 0;
    totalSymbols = 0;
    matchingRatioVec = zeros(1,length(EbN0_dB));

    parfor iter = 1:numIterations
        inv_dset = dSet.^(-1);
        % Generate Random Symbols
        symbols = randi([0 M-1], L, 1); % Random symbols
        % Modulate Symbols (M-QAM)
        modulatedSignal = qammod(symbols, M, 'UnitAveragePower', true);
        
        % Placeholder for Channel (configurable by user)
        delay_p = unifrnd(0,tau_max,[1,K]);
        IT = randn(1,K) + 1i*randn(1,K);
        CR = exp(-2i*pi*rand(1,K));
        %CR = exp(-2i*pi*delay_p*delta_f);
        %CR = exp(2i*pi*delay_p*delta_f);
        channel = IT * transpose(CR).^([0:L-1]); 
        
        % Add AWGN
        % noise = sqrt(noisePower) * (randn(L, 1) + 1i * randn(L, 1));
        receivedSignal = awgn(channel .* modulatedSignal,EbN0_dB(idx),'measured');
        
        % No Channel Estimation (Equalization using true channel)
        equalizedSignal_noEstimation = receivedSignal ./ channel;
        demodulatedSymbols_noEstimation = qamdemod(equalizedSignal_noEstimation, M, 'UnitAveragePower', true);
        numErrors_noEstimation = numErrors_noEstimation + sum(symbols ~= demodulatedSymbols_noEstimation);
        
        % With Channel Estimation (Assume imperfect estimation)
        estimatedChannel1 = awgn(channel,10,"measured"); % 10% error
        % epsilon1 = 0.10; % 10% error
        % noise_std1 = epsilon1 * norm(channel) / sqrt(numel(channel));
        % noise1 = noise_std1 * (randn(size(channel)) + 1i*randn(size(channel))) / sqrt(2);
        % estimatedChannel1 = channel + noise1;
        equalizedSignal_withEstimation1 = receivedSignal ./ estimatedChannel1;
        demodulatedSymbols_withEstimation1 = qamdemod(equalizedSignal_withEstimation1, M, 'UnitAveragePower', true);
        numErrors_withEstimation1 = numErrors_withEstimation1 + sum(symbols ~= demodulatedSymbols_withEstimation1);
        
        % estimatedChannel2 = awgn(channel,20,"measured"); % 1% error
        % % epsilon2 = 0.01; % 1% error
        % % noise_std2 = epsilon2 * norm(channel) / sqrt(numel(channel));
        % % noise2 = noise_std2 * (randn(size(channel)) + 1i*randn(size(channel))) / sqrt(2);
        % % estimatedChannel2 = channel + noise2;
        % equalizedSignal_withEstimation2 = receivedSignal ./ estimatedChannel2;
        % demodulatedSymbols_withEstimation2 = qamdemod(equalizedSignal_withEstimation2, M, 'UnitAveragePower', true);
        % numErrors_withEstimation2 = numErrors_withEstimation2 + sum(symbols ~= demodulatedSymbols_withEstimation2);
        
        estimatedChannel3 = awgn(channel,8,"measured"); % 15% error
        % epsilon3 = 0.15; % 15% error
        % noise_std3 = epsilon3 * norm(channel) / sqrt(numel(channel));
        % noise3 = noise_std3 * (randn(size(channel)) + 1i*randn(size(channel))) / sqrt(2);
        % estimatedChannel3 = channel + noise3;
        equalizedSignal_withEstimation3 = receivedSignal ./ estimatedChannel3;
        demodulatedSymbols_withEstimation3 = qamdemod(equalizedSignal_withEstimation3, M, 'UnitAveragePower', true);
        numErrors_withEstimation3 = numErrors_withEstimation3 + sum(symbols ~= demodulatedSymbols_withEstimation3);
        
        estimatedChannel4 = awgn(channel,7.5,"measured"); % 20% error
        % epsilon4 = 0.20; % 20% error
        % noise_std4 = epsilon4 * norm(channel) / sqrt(numel(channel));
        % noise4 = noise_std4 * (randn(size(channel)) + 1i*randn(size(channel))) / sqrt(2);
        % estimatedChannel4 = channel + noise4;
        equalizedSignal_withEstimation4 = receivedSignal ./ estimatedChannel4;
        demodulatedSymbols_withEstimation4 = qamdemod(equalizedSignal_withEstimation4, M, 'UnitAveragePower', true);
        numErrors_withEstimation4 = numErrors_withEstimation4 + sum(symbols ~= demodulatedSymbols_withEstimation4);
        
        totalSymbols = totalSymbols + L;


        % proposed
        symbols = symbols +1;
        d = dSet(symbols);
        d=d/norm(d);
        %IT = randn(1,K) + 1i*randn(1,K);
        %CR = exp(2i*pi*rand(1,K));
        %CR = exp(-2i*pi*delay_p*delta_f);
        %h = IT * transpose(CR).^([0:L-1]);
        h=channel;
        y = awgn(h.*d, EbN0_dB(idx), "measured");
        [~,~,~,H] = makeHankel(h); % multi-path channel
        [~,~,~,D] = makeHankel(d); % data symbols
        [U,S,V,Y] = makeHankel(y); % received signal
        
        optD = D.^(-1); % optimal decoder
        indexCombinations = [];
        [indexCombinations{1:L}] = ndgrid(1:M);  % Create grid of all combinations
        indexCombinations = reshape(cat(L+1, indexCombinations{:}), [], L);  
        % Reshape into combination list (M^L)
        % indexCombinations = fliplr(indexCombinations);
        numComb = size(indexCombinations, 1);
        detLevel1 = zeros(numComb,1);

        for kk = 1:numComb
            selmatrixA = zeros(M, L);
    
            for col = 1:L
                rowIdx = indexCombinations(kk, col);  % Get the row index for the current combination
                selmatrixA(rowIdx, col) = 1;
            end
    
            search_d = inv_dset*selmatrixA; 
            [U,S,V,search_D] = makeHankel(search_d); % search_D=X_hat
            H_hat = search_D.*Y; % H_hat
            %detLevel1(kk) = abs(det(H_hat*H_hat'));
            [U1,S1,V1]=svd(H_hat);
            s=diag(S1);
            rank1_H= U1(:,1:1)*S1(1:1,1:1)*V1(:,1:1)';
            detLevel1(kk)=norm(H_hat-rank1_H,'fro')^2;
            
        end

        [ss, ii] = min(detLevel1); % find the index of minimum values (minimum determinant index)
        array1 = convertToBinaryVector(indexCombinations(ii,:)-1,M); % 
        array2 = convertToBinaryVector(symbols-1,M); % 
        matchingCountMatrix(iter) = sum(array1 == array2);
    end

    % Compute BER
    matchingRatioVec = mean(matchingCountMatrix) / (L*log2(M));
    BER_proposed(idx) = 1 - matchingRatioVec;
    
    BER_noEstimation(idx) = sum(numErrors_noEstimation) / (log2(M)*totalSymbols);
    BER_withEstimation1(idx) = sum(numErrors_withEstimation1) / (log2(M)*totalSymbols);
    %BER_withEstimation2(idx) = sum(numErrors_withEstimation2) / (log2(M)*totalSymbols);
    BER_withEstimation3(idx) = sum(numErrors_withEstimation3) / (log2(M)*totalSymbols);
    BER_withEstimation4(idx) = sum(numErrors_withEstimation4) / (log2(M)*totalSymbols);
end
BER_proposed       = min(BER_proposed,1);
BER_noEstimation   = min(BER_noEstimation,1);
BER_withEstimation1 = min(BER_withEstimation1,1);
%BER_withEstimation2 = min(BER_withEstimation2,1);
BER_withEstimation3 = min(BER_withEstimation3,1);
BER_withEstimation4 = min(BER_withEstimation4,1);
%% Plot Results
figure(1)
semilogy(EbN0_dB, BER_proposed, 'b-square');
hold on;
semilogy(EbN0_dB, BER_noEstimation, 'r-o');
%semilogy(EbN0_dB, BER_withEstimation2, 'g-^');
semilogy(EbN0_dB, BER_withEstimation1, 'c-s');
semilogy(EbN0_dB, BER_withEstimation3, 'm-+');
semilogy(EbN0_dB, BER_withEstimation4, 'k-d');
xlabel('E_b/N_0 [dB]')
ylabel('BER')
xlim([0, 25]);
legend("Proposed method","Perfect CSIR",'Imperfect CSIR (10% error)','Imperfect CSIR (15% error)','Imperfect CSIR (20% error)','Location','best')
grid on;
ylim([1e-5,1])
