% === 1. Open fig file  ===
fig = openfig('C:\Users\kimyj\OneDrive - 동국대학교\바탕 화면\HankelDecoder_code\BER_a.fig');
ax = gca;
lines = findall(ax, 'Type', 'Line');  

BER_k3 = lines(1).YData;  a_k3 = lines(1).XData;
BER_k2 = lines(2).YData;  a_k2 = lines(2).XData;
BER_k1 = lines(3).YData;  a_k1 = lines(3).XData;

transform_a = @(a) tanh(3 * (a - 1)); 

a_min = min(a_k1);  
a_max = max(a_k1);  
delta = (a_max - a_min) / 2;


x = linspace(-1, 1, 50);  
a_new = 1 + sign(x) .* delta .* (exp(abs(x)) - 1) / (exp(1) - 1);  

BER_k1_new = interp1(a_k1, BER_k1, a_new, 'pchip');
BER_k2_new = interp1(a_k2, BER_k2, a_new, 'pchip');
BER_k3_new = interp1(a_k3, BER_k3, a_new, 'pchip');

figure;
semilogy(a_new, BER_k1_new, '-s', 'LineWidth', 2); hold on;
semilogy(a_new, BER_k2_new, '-o', 'LineWidth', 2);
semilogy(a_new, BER_k3_new, ':^', 'LineWidth', 2);

xlabel(['$|C|$'], ...
        'Interpreter', 'latex', 'FontSize', 11);



ylabel('BER'); grid on;
xlim([min(a_new),max(a_new)])
legend('$K=1$', '$K=2$', '$K=3$', 'Location', 'best','Interpreter', 'latex');
%title('BER vs. $a$ (nonuniform exp spacing)', 'Interpreter', 'latex');
