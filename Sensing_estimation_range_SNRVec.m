clear; clc; close all;
%% Setting the Proposed Method (Hankel Decoder)

phi = 1/13;%phi(41);

C = exp(1i*2*pi*(phi));
%dSet = C.^[0 1 4 6];
dSet = C.^[0 1 3 9];        
%dSet = dSet / sqrt(mean(abs(dSet).^2));
L=2048;                      
delta_f=60e3; % Bandwidth = delta_f * L
c=3e8; % m/s 
delay_per_bin = 1 / (L * delta_f);       % seconds
range_per_bin = c * delay_per_bin / 2;   % meters

% numbers=primes(L-1-1);
% Rset =  numbers(gcd(numbers,L-1)==1);

%% 
% Configuration Parameters
M = 4; % Modulation order 
SNR_dB = -20:2:-4; 

%MSE_per_snr = zeros(size(EbN0_dB));
fail_rate = zeros(size(SNR_dB));
numIterations = 1e3; % Number of iterations

R_real_snr = [];
R_est_snr = [];

T=1; % target

K = 1; % multipaths

Kfactor=6;
tau_max = 200e-9;
R_max=40;

RMSE_per_snr = zeros(size(SNR_dB));
RMSE_zc = zeros(size(SNR_dB));
RMSE_m = zeros(size(SNR_dB));



% Simulation Loop over Eb/N0
for idx = 1:length(SNR_dB)
    rng(idx);
    %EbN0 = 10^(EbN0_dB(idx) / 10); % Linear Eb/N0
    R_real = zeros(numIterations,1);
    R_est = zeros(numIterations,1);
    R_hat_zc = zeros(numIterations,1);
    R_hat_m = zeros(numIterations,1);
    

    parfor iter = 1:numIterations

        H_hat=[];
        
        
        % Generate Random Symbols
        symbols = randi([0 M-1], L, 1); % Random symbols
        
        %R=R_max*rand(1);
        R=unifrnd(0,R_max,[1,T]);
        
        delay_los=(2*R)/c;
   
        IT = randn(1,K-1) + 1i*randn(1,K-1);
        delay_nlos =  rand(1,K-1)*tau_max;
        CR = exp(-2i*pi*delta_f*delay_nlos);       
        channel_NLoS = IT * transpose(CR).^(0:L-1);

        %delay_los=(2*R)./c;
        IT_d = randn(1,1) + 1i*randn(1,1);
        CR_d = exp(-2i*pi*delta_f*delay_los);
        channel_d = IT_d * transpose(CR_d).^(0:L-1);
        channel_d = sqrt(10^(Kfactor/10))*channel_d;


        h=channel_NLoS+channel_d;

        %% proposed
        symbols = symbols +1;
        d = dSet(symbols);
       
        y = awgn(h.*d, SNR_dB(idx), "measured");
        

        y_real = ifft(y); 
        d_ifft = ifft(d); % real transmitted signal

        
        [cc_val, cc_lags] = xcorr(y_real,d_ifft);
        [~, max_idx] = max(abs(cc_val));

        lag_est = cc_lags(max_idx); 
        
        R_hat = lag_est * range_per_bin;

        R_real(iter)= (R);
        R_est(iter)= (R_hat);

%% ZC sequence
        %r = Rset(randi(length(Rset),1));
        L2=L-1;
        % h2=h(1:L2);

        IT = randn(1,K-1) + 1i*randn(1,K-1);
        delay_nlos =  rand(1,K-1)*tau_max;
        CR = exp(-2i*pi*delta_f*delay_nlos);       
        channel_NLoS = IT * transpose(CR).^(0:L2-1);

        %delay_los=(2*R)./c;
        IT_d = randn(1,1) + 1i*randn(1,1);
        CR_d = exp(-2i*pi*delta_f*delay_los);
        channel_d = IT_d * transpose(CR_d).^(0:L2-1);
        channel_d = sqrt(10^(Kfactor/10))*channel_d;


        h2=channel_NLoS+channel_d;
        r=5;
        xx = zadoffChuSeq(r,L2);
        %xx= xx / sqrt(mean(abs(xx).^2));
        y2 = awgn(h2.*(xx.'), SNR_dB(idx), "measured");

        xx_ifft=ifft(xx);
        y2_real=ifft(y2);

        [cc_val2, cc_lags2] = xcorr(y2_real,xx_ifft);
        [~, max_idx2] = max(abs(cc_val2));

        lag_est2 = cc_lags2(max_idx2); 
        
        R_hat2 = lag_est2 * range_per_bin;

        R_hat_zc(iter) = (R_hat2);


        %% M-seq.
        m = mlseq(L2);
        %m= m / sqrt(mean(abs(m).^2));
        %m=transpose(m);
        y3 = awgn(h2.*(m.'),SNR_dB(idx),'measured');

        m_ifft=ifft(m);
        y3_real=ifft(y3);

        [cc_val3, cc_lags3] = xcorr(y3_real,m_ifft);
        [~, max_idx3] = max(abs(cc_val3));

        lag_est3 = cc_lags3(max_idx3); 
        
        R_hat3 = lag_est3 * range_per_bin;

        R_hat_m(iter) = (R_hat3);

        
        
    end
 
RMSE_per_snr(idx) = rmse(R_est,R_real);
RMSE_zc(idx) = rmse(R_hat_zc,R_real);  
RMSE_m(idx) = rmse(R_hat_m,R_real);


end


%% Plot Results
figure(1);
semilogy(SNR_dB, (RMSE_per_snr)./(R_max), '-s', 'LineWidth', 2);
hold on;
semilogy(SNR_dB, (RMSE_zc)./(R_max), '--o', 'LineWidth', 2);
hold on
semilogy(SNR_dB, (RMSE_m)./(R_max), '-+', 'LineWidth', 2);
hold on;
xlabel('SNR [dB]');
xticks(SNR_dB)
ylabel('Normalized range error (R_{rmse}/R_{max})');
legend('Proposed method (DFS-GR)', 'Baseline 1 (ZC-sequence)', 'Baseline 2 (m-sequence)', 'Location', 'best', 'fontsize',10);
%ylim([1e-2,10])
grid on;
%semilogy(SNR_dB, sqrt(CRB)./R_max)
% delta_R = c / (2 * (L) * delta_f);       % Resolution in meters
% normalized_resolution = delta_R / R_max;  % Normalized resolution
% 
% yline(normalized_resolution, 'k--', 'LineWidth', 2);

%%
% figure(2);
% semilogy(EbN0_dB, sqrt(MSE_per_snr), '-s', 'LineWidth', 2);
% hold on;
% semilogy(EbN0_dB, sqrt(MSE_zc), '--o', 'LineWidth', 2);
% hold on
% semilogy(EbN0_dB, sqrt(MSE_m), '-+', 'LineWidth', 2);
% hold on;
% xlabel('SNR [dB]');
% xticks(EbN0_dB)
% ylabel('Range error (R_{rmse})');
% legend('Proposed (Golomb ruler-based sequence)', 'Baseline 1 (ZC-sequence)', 'Baseline 2 (m-sequence)', 'Location', 'best', 'fontsize',10);
% %ylim([1e-2,10])
% grid on;
