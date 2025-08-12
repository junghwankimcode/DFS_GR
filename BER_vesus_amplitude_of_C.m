clear; clc;
%% Setting the Proposed Method (Hankel Decoder)

mvec=[0.7:0.05:1.3];
%phi = 0.01:0.01:1;
phi = 1/13;
%C = m*exp(1i*2*pi*(phi));
C = exp(1i*2*pi*(phi));
%dSet = C.^[0 1 4 6];
dSet = C.^[0 1 3 9];
delta_f = 15*1e3;  
tau_max = 1e-6;


%% M-QAM BER Simulator with Channel Estimation
% Configuration Parameters
M = 4; % Modulation order 
%EbN0_dB = 0:5:25; % Range of Eb/N0 values in dB
EbN0 = 20;
numIterations = 2e3; % Number of iterations

L = 8; % Length of received signal per iteration
%K = 1;
Kvec=[1,2,3];


BER_proposed_K=zeros(length(Kvec),length(mvec));



% Simulation Loop over Eb/N0
for k=1:length(Kvec)
    K=Kvec(k)


for idx = 1:length(mvec)
    m=mvec(idx);
    
    C = m*exp(1i*2*pi*(phi));
    dSet = C.^[0 1 3 9];
    %dSet = dSet / sqrt(mean(abs(dSet).^2));
    %EbN0 = 10^(EbN0_dB(idx) / 10); % Linear Eb/N0
    noisePower = 1 / (2 * log2(M) * EbN0); % Noise power per dimension
    numErrors_noEstimation = 0;
    numErrors_withEstimation1 = 0;
    numErrors_withEstimation2 = 0;
    numErrors_withEstimation3 = 0;
    numErrors_withEstimation4 = 0;
    totalSymbols = 0;
    matchingRatioVec = zeros(1,length(mvec));
    %BER_proposed = zeros(1,length(mvec));
    parfor iter = 1:numIterations
        inv_dset = dSet.^(-1);
        % Generate Random Symbols
        symbols = randi([0 M-1], L, 1); % Random symbols
        % Modulate Symbols (M-QAM)
        %modulatedSignal = qammod(symbols, M, 'UnitAveragePower', true);
        
        % Placeholder for Channel (configurable by user)
        delay_p = unifrnd(0,tau_max,[1,K]);
        IT = randn(1,K) + 1i*randn(1,K);
        %CR = exp(-2i*pi*rand(1,K));
        CR = exp(-2i*pi*delay_p*delta_f);
        %CR = exp(2i*pi*delay_p*delta_f);
        channel = IT * transpose(CR).^([0:L-1]); 
        
        
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
        y = awgn(h.*d, EbN0, "measured");
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
    
    
end

BER_proposed_K(k,:)=BER_proposed


end

%% Plot Results
figure(1)
semilogy(mvec, BER_proposed_K(1,:), '-square','LineWidth',2);
hold on;
semilogy(mvec, BER_proposed_K(2,:), '--*','LineWidth',2);
hold on;
semilogy(mvec, BER_proposed_K(3,:), ':+','LineWidth',2);
hold on;
xlabel('$a$', 'Interpreter', 'latex')
ylabel('BER')
xticks(mvec)
xlim([min(mvec), max(mvec)]);
legend('$K$=1','$K$=2','$K$=3','Location','best','Interpreter', 'latex','fontsize',12)
grid on;
%ylim([1e-5,1])
