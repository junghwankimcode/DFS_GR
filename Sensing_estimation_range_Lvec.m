clear; clc; close all;
%% Setting the Proposed Method (Hankel Decoder)

phi = 1/13;%phi(41);
C = exp(1i*2*pi*(phi));
%dSet = C.^[0 1 4 6];
dSet = C.^[0 1 3 9];                            
delta_f=60e3; % Bandwidth = delta_f * L
c=3e8; % m/s 
% delay_per_bin = 1 / (L * delta_f);       % seconds
% range_per_bin = c * delay_per_bin / 2;   % meters

% numbers=primes(L-1-1);
% Rset =  numbers(gcd(numbers,L-1)==1);

%% 
% Configuration Parameters
M = 4; % Modulation order 
LVec = [128, 256, 512, 1024, 2048]; 
SNR =-10;
%MSE_per_snr = zeros(size(EbN0_dB));
%fail_rate = zeros(size(LVec));
numIterations = 1e3; % Number of iterations

R_real_snr = [];
R_est_snr = [];

T=1; % target

K = 1; % multipaths

Kfactor=6;
tau_max = 200e-9;
R_max=40;

RMSE_per_snr = zeros(size(LVec));
RMSE_zc = zeros(size(LVec));
RMSE_m = zeros(size(LVec));


CRLB=zeros(size(LVec));

% Simulation Loop over Eb/N0
for idx = 1:length(LVec)
    rng(idx);
    %EbN0 = 10^(EbN0_dB(idx) / 10); % Linear Eb/N0
    R_real = zeros(numIterations,1);
    R_est = zeros(numIterations,1);
    R_hat_zc = zeros(numIterations,1);
    R_hat_m = zeros(numIterations,1);
    
    L=LVec(idx);
    delay_per_bin = 1 / (L * delta_f);       % seconds
    range_per_bin = c * delay_per_bin / 2;   % meters
    delay_per_bin2 = 1 / ((L-1) * delta_f);       % seconds
    range_per_bin2 = c * delay_per_bin2 / 2;   % meters

    parfor iter = 1:numIterations

        H_hat=[];
        
        
        % Generate Random Symbols
        symbols = randi([0 M-1], L, 1); % Random symbols
        
        R=R_max*rand(1);%unifrnd(0,R_max,[1,T]);
        
        delay_los=(2*R)/c;
   
        IT = randn(1,K-1) + 1i*randn(1,K-1);
        delay_nlos = rand(1,K-1)*tau_max;
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
       
        y = awgn(h.*d, SNR, "measured");
        

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
        delay_nlos = rand(1,K-1)*tau_max;
        CR = exp(-2i*pi*delta_f*delay_nlos);       
        channel_NLoS = IT * transpose(CR).^(0:L2-1);

        %delay_los=(2*R)./c;
        IT_d = randn(1,1) + 1i*randn(1,1);
        CR_d = exp(-2i*pi*delta_f*delay_los);
        channel_d = IT_d * transpose(CR_d).^(0:L2-1);
        channel_d = sqrt(10^(Kfactor/10))*channel_d;


        h2=channel_NLoS+channel_d;
        numbers=primes(L2-1);
        r_cand= numbers(gcd(numbers,L2)==1);
        r = r_cand(randi(length(r_cand)));
        xx = zadoffChuSeq(r,L2);
        %xx= xx / sqrt(mean(abs(xx).^2));
        y2 = awgn(h2.*(xx.'), SNR, "measured");

        xx_ifft=ifft(xx);
        y2_real=ifft(y2);

        [cc_val2, cc_lags2] = xcorr(y2_real,xx_ifft);
        [~, max_idx2] = max(abs(cc_val2));

        lag_est2 = cc_lags2(max_idx2); 
        
        R_hat2 = lag_est2 * range_per_bin2;

        R_hat_zc(iter) = (R_hat2);


        %% M-seq.
        m = mlseq(L2);
        %m= m / sqrt(mean(abs(m).^2));
        %m=transpose(m);
        y3 = awgn(h2.*(m.'), SNR,'measured');

        m_ifft=ifft(m);
        y3_real=ifft(y3);

        [cc_val3, cc_lags3] = xcorr(y3_real,m_ifft);
        [~, max_idx3] = max(abs(cc_val3));

        lag_est3 = cc_lags3(max_idx3); 
        
        R_hat3 = lag_est3 * range_per_bin2;

        R_hat_m(iter) = (R_hat3);

        
        
    end
 
RMSE_per_snr(idx) = rmse(R_est,R_real);
RMSE_zc(idx) = rmse(R_hat_zc,R_real);  
RMSE_m(idx) = rmse(R_hat_m,R_real);


end


%% Plot Results
% figure(1);
% semilogy(LVec, (RMSE_per_snr)./(R_max), '-s', 'LineWidth', 2);
% hold on;
% semilogy(LVec, (RMSE_zc)./(R_max), '--o', 'LineWidth', 2);
% hold on
% semilogy(LVec, (RMSE_m)./(R_max), '-+', 'LineWidth', 2);
% hold on;
% xlabel('$L$','Interpreter', 'latex');
% xticks(LVec)
% ylabel('Normalized range error (R_{rmse}/R_{max})');
% legend('Proposed (Golomb ruler-based sequence)', 'Baseline 1 (ZC-sequence (Length = $L$-1))', 'Baseline 2 (m-sequence (Length = $L$-1))', ...
%     'Location', 'best', 'fontsize',10, 'Interpreter','latex');
% %ylim([1e-2,10])
% xlim([min(LVec),max(LVec)])
% grid on;

%%
figure(2);
clf;
x_idx = 1:length(LVec); % 등간격 x축 인덱스

semilogy(x_idx, RMSE_per_snr./R_max, '-s', 'LineWidth', 2); hold on;
semilogy(x_idx, RMSE_zc./R_max, '--o', 'LineWidth', 2); hold on;
semilogy(x_idx, RMSE_m./R_max, '-+', 'LineWidth', 2); hold on;

xticks(x_idx);
xticklabels(LVec);  % 진짜 L 값들을 레이블로 표시
xlabel('$L$', 'Interpreter', 'latex');
ylabel('Normalized range error (R_{rmse}/R_{max})');
legend({'Proposed (Golomb ruler-based sequence)', ...
        'Baseline 1 (ZC-sequence)', ...
        'Baseline 2 (m-sequence)'}, ...
        'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
grid on;
ylim([1e-2,10])