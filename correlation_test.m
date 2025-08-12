clear; clc; close all


phi = 1/13;%phi(41);
%C = m*exp(1i*2*pi*(phi));
C = exp(1i*2*pi*(phi));
%dSet = C.^[0 1 4 6];
dSet = C.^[0 1 3 9];  
dSet = dSet / sqrt(mean(abs(dSet).^2));
L=121;                      
M=4;
numbers=primes(L-1);
Rset =  numbers(gcd(numbers,L)==1);

symbols = randi([0 M-1], L, 1); % Random symbols
symbols = symbols +1;
d = dSet(symbols);

%R = Rset(randi(length(Rset),1));
R1=5;
R2=7;
R3=113;
xx = zadoffChuSeq(R1,L);
xx=xx/sqrt(mean(abs(xx).^2)^2);
xx2 = zadoffChuSeq(R2,L);
xx2=xx2/sqrt(mean(abs(xx2).^2)^2);
xx3 = zadoffChuSeq(R3,L);
xx3=xx3/sqrt(mean(abs(xx3).^2)^2);
%SNR = 10;

%d = awgn(d,SNR,'measured');
%xx = awgn(xx,SNR,'measured');

%d=d/norm(d);
%xx=xx/norm(xx);

d1=awgn(d,0,"measured");
x1=awgn(xx,0,"measured");
x2=awgn(xx2,0,"measured");
x3=awgn(xx3,0,'measured');

[ab1,lag1]=xcorr(d1,d1);
[ab2,lag2]=xcorr(x1,x1);
[ab3,lag3]=xcorr(x2,x2);
[ab4,lag4]=xcorr(x3,x3);

figure(1)
subplot(2,1,1)
stem(lag1,abs(ab1)./ max(abs(ab1)),'filled')
hold on;
stem(lag2,abs(ab2)./ max(abs(ab2)))
grid on;
stem(lag3,abs(ab3)./ max(abs(ab3)))
grid on;
stem(lag3,abs(ab4)./ max(abs(ab4)))
grid on;
xlabel('lag')
ylabel('Magnitude')
%ylim([0,1])
%legend('Golomb ruler-based sequence','ZC-sequence','Location','best','fontsize',10)
legend('Proposed method (DFS-GR)','ZC-sequence (root index=5)','ZC-sequence (root index=7)','ZC-sequence (root index=113)', ...
    'Location','best','fontsize',10)

% golomb_pslr=computePSLR(d)
% zc_pslr1=computePSLR((xx))


%%
%  clear; clc;
% 
% m = 0.5:0.01:2;
% m = m(2);
% phi = 0.01:0.01:1;
% phi = 1/13;%phi(41);
% %C = m*exp(1i*2*pi*(phi));
% C = exp(1i*2*pi*(phi));
% %dSet = C.^[0 1 4 6];
% dSet = C.^[0 1 3 9];       
% dSet = dSet / sqrt(mean(abs(dSet).^2));
% L=121;                      
% M=4;
% numbers=primes(L-1);
% Rset =  numbers(gcd(numbers,L)==1);
% 
% symbols = randi([0 M-1], L, 1); % Random symbols
% symbols = symbols +1;
% d = dSet(symbols);
% 
% % R = Rset(randi(length(Rset),1));
% R1=2;
% R2=3;
% R3=5;
% R4=7;
% R5=113;
% xx1 = zadoffChuSeq(R1,L);
% xx1=xx1/sqrt(mean(abs(xx1).^2)^2);
% xx2 = zadoffChuSeq(R2,L);
% xx2=xx2/sqrt(mean(abs(xx2).^2)^2);
% xx3 = zadoffChuSeq(R3,L);
% xx3=xx3/sqrt(mean(abs(xx3).^2)^2);
% xx4 = zadoffChuSeq(R4,L);
% xx4=xx4/sqrt(mean(abs(xx4).^2)^2);
% xx5 = zadoffChuSeq(R5,L);
% xx5=xx5/sqrt(mean(abs(xx5).^2)^2);
% 
% SNRVec=[0:1:30];
% 
% % d=d/norm(d);
% % xx1=xx1/norm(xx1);
% % xx2=xx2/norm(xx2);
% 
% 
% golomb_pslr_vec_f=[];
% zc_pslr_vec_f1=[];
% zc_pslr_vec_f2=[];
% zc_pslr_vec_f3=[];
% zc_pslr_vec_f4=[];
% zc_pslr_vec_f5=[];
% 
% for i = 1:length(SNRVec)
% 
%     golomb_pslr_vec=[];
%     zc_pslr_vec1=[];
%     zc_pslr_vec2=[];
%     zc_pslr_vec3=[];
%     zc_pslr_vec4=[];
%     zc_pslr_vec5=[];
% 
%     for k=1:1e3
% 
%         % symbols = randi([0 M-1], L, 1); % Random symbols
%         % symbols = symbols +1;
%         % d = dSet(symbols);
%         % 
%         % %R = Rset(randi(length(Rset),1));
%         % R=5;
%         % xx = zadoffChuSeq(R,L);
% 
%         SNR=SNRVec(i);
% 
%         golomb_pslr=0;
%         zc_pslr1=0;
%         zc_pslr2=0;zc_pslr3=0;zc_pslr4=0;zc_pslr5=0;
% 
%         d_noisy = awgn(d,SNR,'measured');
%         xx_noisy1 = awgn(xx1,SNR,'measured');
%         xx_noisy2 = awgn(xx2,SNR,'measured');
%         xx_noisy3 = awgn(xx3,SNR,'measured');
%         xx_noisy4 = awgn(xx4,SNR,'measured');
%         xx_noisy5 = awgn(xx5,SNR,'measured');
% 
% 
% 
%         golomb_pslr=computePSLR(d_noisy);
%         zc_pslr1=computePSLR((xx_noisy1));
%         zc_pslr2=computePSLR((xx_noisy2));
%         zc_pslr3=computePSLR((xx_noisy3));
%         zc_pslr4=computePSLR((xx_noisy4));
%         zc_pslr5=computePSLR((xx_noisy5));
% 
%         golomb_pslr_vec=[golomb_pslr_vec, golomb_pslr];
%         zc_pslr_vec1=[zc_pslr_vec1, zc_pslr1];
%         zc_pslr_vec2=[zc_pslr_vec2, zc_pslr2];
%         zc_pslr_vec3=[zc_pslr_vec3, zc_pslr3];
%         zc_pslr_vec4=[zc_pslr_vec4, zc_pslr4];
%         zc_pslr_vec5=[zc_pslr_vec5, zc_pslr5];
%     end
% 
%     golomb_pslr_vec_f=[golomb_pslr_vec_f mean(golomb_pslr_vec)];
%     zc_pslr_vec_f1=[zc_pslr_vec_f1 mean(zc_pslr_vec1)];
%     zc_pslr_vec_f2=[zc_pslr_vec_f2 mean(zc_pslr_vec2)];
%     zc_pslr_vec_f3=[zc_pslr_vec_f3 mean(zc_pslr_vec3)];
%     zc_pslr_vec_f4=[zc_pslr_vec_f4 mean(zc_pslr_vec4)];
%     zc_pslr_vec_f5=[zc_pslr_vec_f5 mean(zc_pslr_vec5)];
% 
% end
% 
% figure(2)
% % subplot(2,1,2)
% plot(SNRVec, golomb_pslr_vec_f,'-bs', 'LineWidth', 2)
% hold on
% % plot(SNRVec, zc_pslr_vec_f1,'-o', 'LineWidth', 2)
% % hold on
% % plot(SNRVec, zc_pslr_vec_f2,'-+', 'LineWidth', 2)
% % hold on
% plot(SNRVec, zc_pslr_vec_f3,'-d', 'LineWidth', 2)
% hold on
% plot(SNRVec, zc_pslr_vec_f4,'-^', 'LineWidth', 2)
% hold on
% plot(SNRVec, zc_pslr_vec_f5,'->', 'LineWidth', 2)
% hold on
% grid on
% xlabel('SNR [dB]')
% ylabel('Peak-to-sidelobe ratio [dB]')
% xticks([0:5:30])
% legend('Golomb ruler-based sequence','ZC-sequence (root index=5)','ZC-sequence (root index=7)','ZC-sequence (root index=113)', ...
%     'Location','best','fontsize',10)


%%
%clear; clc; close all

m = 0.5:0.01:2;
m = m(67);
phi = 0.01:0.01:1;
phi = 1/13;%phi(41);
%C = m*exp(1i*2*pi*(phi));
C = exp(1i*2*pi*(phi));
%dSet = C.^[0 1 4 6];
dSet = C.^[0 1 3 9];  
dSet = dSet / sqrt(mean(abs(dSet).^2));
L=121;                      
M=4;
numbers=primes(L-1);
Rset =  numbers(gcd(numbers,L)==1);

symbols = randi([0 M-1], L, 1); % Random symbols
symbols = symbols +1;
d = dSet(symbols);

%R = Rset(randi(length(Rset),1));
R1=5;
R2=7;
R3=113;
xx = zadoffChuSeq(R1,L);
xx=xx/sqrt(mean(abs(xx).^2)^2);
xx2 = zadoffChuSeq(R2,L);
xx2=xx2/sqrt(mean(abs(xx2).^2)^2);
xx3 = zadoffChuSeq(R3,L);
xx3=xx3/sqrt(mean(abs(xx3).^2)^2);
%SNR = 10;

%d = awgn(d,SNR,'measured');
%xx = awgn(xx,SNR,'measured');

%d=d/norm(d);
%xx=xx/norm(xx);

d1=awgn(d,0,"measured");
x1=awgn(xx,0,'measured');
x2=awgn(xx2,0,"measured");
x3=awgn(xx3,0,"measured");

% [ab1,lag1]=xcorr(d1,d1);
% [ab2,lag2]=xcorr(x1,x1);
% [ab3,lag3]=xcorr(x2,x2);
% [ab4,lag4]=xcorr(x3,x3);
ab1=circular_autocorr(d1);
ab2=circular_autocorr(x1);
ab3=circular_autocorr(x2);
ab4=circular_autocorr(x3);


x=0:L-1;
figure(1)
subplot(2,1,2)
stem(x,abs(ab1)./ max(abs(ab1)),'filled','LineWidth',2)
hold on;
stem(x,abs(ab2)./ max(abs(ab2)),'LineWidth',2)
grid on;
stem(x,abs(ab3)./ max(abs(ab3)),'LineWidth',2)
grid on;
stem(x,abs(ab4)./ max(abs(ab4)),'LineWidth',2)
grid on;
xlabel('Circular shift')
ylabel('Magnitude')
%ylim([0,1])
xlim([0 L-1])
%legend('Golomb ruler-based sequence','ZC-sequence','Location','best','fontsize',10)
legend('Proposed method (DFS-GR)','ZC-sequence (root index=5)','ZC-sequence (root index=7)','ZC-sequence (root index=113)', ...
    'Location','best','fontsize',10)


%%
clear; clc;

m = 0.5:0.01:2;
m = m(2);
phi = 0.01:0.01:1;
phi = 1/13;%phi(41);
%C = m*exp(1i*2*pi*(phi));
C = exp(1i*2*pi*(phi));
%dSet = C.^[0 1 4 6];
dSet = C.^[0 1 3 9];       
%dSet = dSet / sqrt(mean(abs(dSet).^2));
L=121;                      
M=4;
numbers=primes(L-1);
Rset =  numbers(gcd(numbers,L)==1);

% symbols = randi([0 M-1], L, 1); % Random symbols
% symbols = symbols +1;
% d = dSet(symbols);
% 
% % R = Rset(randi(length(Rset),1));
% R1=2;
% R2=3;
% R3=5;
% R4=7;
% R5=113;
% xx1 = zadoffChuSeq(R1,L);
% %xx1=xx1/sqrt(mean(abs(xx1).^2)^2);
% xx2 = zadoffChuSeq(R2,L);
% %xx2=xx2/sqrt(mean(abs(xx2).^2)^2);
% xx3 = zadoffChuSeq(R3,L);
% %xx3=xx3/sqrt(mean(abs(xx3).^2)^2);
% xx4 = zadoffChuSeq(R4,L);
% %xx4=xx4/sqrt(mean(abs(xx4).^2)^2);
% xx5 = zadoffChuSeq(R5,L);
% %xx5=xx5/sqrt(mean(abs(xx5).^2)^2);

SNRVec=-15:1:15;

% d=d/norm(d);
% xx1=xx1/norm(xx1);
% xx2=xx2/norm(xx2);


golomb_pslr_vec_f=[];
zc_pslr_vec_f1=[];
zc_pslr_vec_f2=[];
zc_pslr_vec_f3=[];
zc_pslr_vec_f4=[];
zc_pslr_vec_f5=[];

golomb_pslr_vec_f_c=[];
zc_pslr_vec_f1_c=[];
zc_pslr_vec_f2_c=[];
zc_pslr_vec_f3_c=[];
zc_pslr_vec_f4_c=[];
zc_pslr_vec_f5_c=[];

for i = 1:length(SNRVec)

    golomb_pslr_vec=[];
    zc_pslr_vec1=[];
    zc_pslr_vec2=[];
    zc_pslr_vec3=[];
    zc_pslr_vec4=[];
    zc_pslr_vec5=[];

    golomb_pslr_vec_c=[];
    zc_pslr_vec1_c=[];
    zc_pslr_vec2_c=[];
    zc_pslr_vec3_c=[];
    zc_pslr_vec4_c=[];
    zc_pslr_vec5_c=[];

    for k=1:1e3

        % C = exp(1i*2*pi*(phi));
        % %dSet = C.^[0 1 4 6];
        % dSet = C.^[0 1 3 9];       
        % dSet = dSet / sqrt(mean(abs(dSet).^2));
        % L=121;                      
        % M=4;
        % numbers=primes(L-1);
        % Rset =  numbers(gcd(numbers,L)==1);
        % 
        symbols = randi([0 M-1], L, 1); % Random symbols
        symbols = symbols +1;
        d = dSet(symbols);

        % R = Rset(randi(length(Rset),1));
        R1=2;
        R2=3;
        R3=5;
        R4=7;
        R5=113;
        xx1 = zadoffChuSeq(R1,L);
        %xx1=xx1/sqrt(mean(abs(xx1).^2)^2);
        xx2 = zadoffChuSeq(R2,L);
        %xx2=xx2/sqrt(mean(abs(xx2).^2)^2);
        xx3 = zadoffChuSeq(R3,L);
        %xx3=xx3/sqrt(mean(abs(xx3).^2)^2);
        xx4 = zadoffChuSeq(R4,L);
        %xx4=xx4/sqrt(mean(abs(xx4).^2)^2);
        xx5 = zadoffChuSeq(R5,L);
        %xx5=xx5/sqrt(mean(abs(xx5).^2)^2);


        SNR=SNRVec(i);

        golomb_pslr=0;
        zc_pslr1=0;
        zc_pslr2=0;zc_pslr3=0;zc_pslr4=0;zc_pslr5=0;

        d_noisy = awgn(d,SNR,'measured');
        xx_noisy1 = awgn(xx1,SNR,'measured');
        xx_noisy2 = awgn(xx2,SNR,'measured');
        xx_noisy3 = awgn(xx3,SNR,'measured');
        xx_noisy4 = awgn(xx4,SNR,'measured');
        xx_noisy5 = awgn(xx5,SNR,'measured');



        golomb_pslr=computePSLR(d_noisy);
        zc_pslr1=computePSLR((xx_noisy1));
        zc_pslr2=computePSLR((xx_noisy2));
        zc_pslr3=computePSLR((xx_noisy3));
        zc_pslr4=computePSLR((xx_noisy4));
        zc_pslr5=computePSLR((xx_noisy5));

        golomb_pslr_c=computePSLR_circular(d_noisy);
        zc_pslr1_c=computePSLR_circular(xx_noisy1);
        zc_pslr2_c=computePSLR_circular(xx_noisy2);
        zc_pslr3_c=computePSLR_circular(xx_noisy3);
        zc_pslr4_c=computePSLR_circular(xx_noisy4);
        zc_pslr5_c=computePSLR_circular(xx_noisy5);

        golomb_pslr_vec=[golomb_pslr_vec, golomb_pslr];
        zc_pslr_vec1=[zc_pslr_vec1, zc_pslr1];
        zc_pslr_vec2=[zc_pslr_vec2, zc_pslr2];
        zc_pslr_vec3=[zc_pslr_vec3, zc_pslr3];
        zc_pslr_vec4=[zc_pslr_vec4, zc_pslr4];
        zc_pslr_vec5=[zc_pslr_vec5, zc_pslr5];

        golomb_pslr_vec_c=[golomb_pslr_vec_c, golomb_pslr_c];
        zc_pslr_vec1_c=[zc_pslr_vec1_c, zc_pslr1_c];
        zc_pslr_vec2_c=[zc_pslr_vec2_c, zc_pslr2_c];
        zc_pslr_vec3_c=[zc_pslr_vec3, zc_pslr3];
        zc_pslr_vec4_c=[zc_pslr_vec4_c, zc_pslr4_c];
        zc_pslr_vec5_c=[zc_pslr_vec5_c, zc_pslr5_c];
    end

    golomb_pslr_vec_f=[golomb_pslr_vec_f mean(golomb_pslr_vec)];
    zc_pslr_vec_f1=[zc_pslr_vec_f1 mean(zc_pslr_vec1)];
    zc_pslr_vec_f2=[zc_pslr_vec_f2 mean(zc_pslr_vec2)];
    zc_pslr_vec_f3=[zc_pslr_vec_f3 mean(zc_pslr_vec3)];
    zc_pslr_vec_f4=[zc_pslr_vec_f4 mean(zc_pslr_vec4)];
    zc_pslr_vec_f5=[zc_pslr_vec_f5 mean(zc_pslr_vec5)];

    golomb_pslr_vec_f_c=[golomb_pslr_vec_f_c mean(golomb_pslr_vec_c)];
    zc_pslr_vec_f1_c=[zc_pslr_vec_f1_c mean(zc_pslr_vec1_c)];
    zc_pslr_vec_f2_c=[zc_pslr_vec_f2_c mean(zc_pslr_vec2_c)];
    zc_pslr_vec_f3_c=[zc_pslr_vec_f3_c mean(zc_pslr_vec3_c)];
    zc_pslr_vec_f4_c=[zc_pslr_vec_f4_c mean(zc_pslr_vec4_c)];
    zc_pslr_vec_f5_c=[zc_pslr_vec_f5_c mean(zc_pslr_vec5_c)];

end

%
figure(2)
subplot(2,1,1)
plot(SNRVec, golomb_pslr_vec_f,'-bs', 'LineWidth', 2)
hold on
% plot(SNRVec, zc_pslr_vec_f1,'-o', 'LineWidth', 2)
% hold on
% plot(SNRVec, zc_pslr_vec_f2,'-+', 'LineWidth', 2)
% hold on
plot(SNRVec, zc_pslr_vec_f3,':d', 'LineWidth', 1)
hold on
plot(SNRVec, zc_pslr_vec_f4,'-.^', 'LineWidth', 1)
hold on
plot(SNRVec, zc_pslr_vec_f5,'-->', 'LineWidth', 1)
hold on
grid on
xlabel('SNR [dB]')
ylabel('Peak-to-sidelobe level ratio [dB]')
xticks([min(SNRVec):5:max(SNRVec)])
legend('Proposed method (DFS-GR)','ZC-sequence (root index=5)','ZC-sequence (root index=7)','ZC-sequence (root index=113)', ...
    'Location','best','fontsize',9,'NumColumns',1)



subplot(2,1,2)
plot(SNRVec, golomb_pslr_vec_f_c,'-bs', 'LineWidth', 2)
hold on
% plot(SNRVec, zc_pslr_vec_f1,'-o', 'LineWidth', 2)
% hold on
% plot(SNRVec, zc_pslr_vec_f2,'-+', 'LineWidth', 2)
% hold on
plot(SNRVec, zc_pslr_vec_f3_c,':d', 'LineWidth', 1)
hold on
plot(SNRVec, zc_pslr_vec_f4_c,'-.^', 'LineWidth', 1)
hold on
plot(SNRVec, zc_pslr_vec_f5_c,'-->', 'LineWidth', 1)
hold on
grid on
xlabel('SNR [dB]')
ylabel('Peak-to-sidelobe level ratio [dB]')
xticks([min(SNRVec):5:max(SNRVec)])
legend('Proposed method (DFS-GR)','ZC-sequence (root index=5)','ZC-sequence (root index=7)','ZC-sequence (root index=113)', ...
    'Location','best','fontsize',9,'NumColumns',1)
