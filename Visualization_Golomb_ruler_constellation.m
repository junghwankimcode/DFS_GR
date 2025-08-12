clear; clc;

% 1) Parameters
phi      = 1/13;                 % base phase
marks    = [0,1,3,9];            % Golomb‐ruler marks
A = [0.95];
C        = A*exp(1i*2*pi*phi);
dSet     = C.^marks;             % your 4‐point constellation
dSet = dSet / sqrt(mean(abs(dSet).^2));
M        = numel(dSet);

L        = 8;                    % symbols per frame
EbN0_dB  = 22;                   % demo SNR
numIter  = 2000;                 % Monte Carlo trials
tau_max  = 1e-6;                % channel spread
delta_f  = 15e3;                % frequency offset

% 2) Build equalized receive cloud (perfect CSIR)
eq_points = [];
for iter = 1:numIter
    idx     = randi([1 M], L, 1);
    d       = dSet(idx);
    delay_p = unifrnd(0, tau_max);
    IT      = randn + 1i*randn;
    h       = IT * exp(-2i*pi*delay_p*delta_f).^(0:L-1);
    rx      = awgn(h.' .* d, EbN0_dB, 'measured');
    eq      = rx ./ (h.');            % perfect equalization
    eq_points = [eq_points; eq(:)];
end

% 3) Compute 6 CW + 6 CCW + full circle ⇒ 13 arcs
angles    = 2*pi*phi * marks;                % absolute phases
pair_list = nchoosek(1:M,2);                 % all 6 i<j pairs
Npairs    = size(pair_list,1);

d_cw  = zeros(Npairs,1);
d_ccw = zeros(Npairs,1);
for k = 1:Npairs
    i = pair_list(k,1);
    j = pair_list(k,2);
    d_cw(k)  = mod( angles(j) - angles(i), 2*pi );
    d_ccw(k) = mod( angles(i) - angles(j), 2*pi );
end
full_cw = 2*pi;
all_arcs  = [d_cw; d_ccw; full_cw];
arc_pairs = [pair_list; flip(pair_list,2); NaN NaN];
arcColors = hsv(numel(all_arcs));

% legend labels for arcs only
arcLabels = cell(numel(all_arcs),1);
for k = 1:numel(all_arcs)
    if k <= Npairs || k <= 2*Npairs
        i = arc_pairs(k,1); j = arc_pairs(k,2);
        arcLabels{k} = sprintf('%d→%d (%.2fπ)', marks(i), marks(j), all_arcs(k)/pi);
    else
        arcLabels{k} = 'full circle (2π)';
    end
end

% 4) Plot
figure;
hold on;

% a) Rx cloud (hidden in legend)
scatter(real(eq_points), imag(eq_points), 8, '.', 'MarkerEdgeAlpha',0.15, 'HandleVisibility','off');

% b) unit‐circle (no legend)
theta = linspace(0,2*pi,360);
plot(cos(theta), sin(theta), 'k--','LineWidth',1,'HandleVisibility','off');

% c) Tx constell + spokes (hidden in legend)
scatter(real(dSet), imag(dSet), 100, 'r+','LineWidth',1, 'HandleVisibility','off');
for i = 1:M
    plot([0 real(dSet(i))], [0 imag(dSet(i))], 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    text(real(dSet(i))*1.15, imag(dSet(i))*1.15, sprintf('%d',marks(i)), 'HorizontalAlignment','center', 'FontSize',12,'Color','g','FontWeight','bold');
end

% d) draw 13 arcs and collect handles
arcHandles = gobjects(numel(all_arcs),1);
for k = 1:numel(all_arcs)
    r    = 2 + 0.2*k;
    i    = arc_pairs(k,1);
    if isnan(i)
        arc_t = theta;
    else
        arc_t = linspace(angles(i), angles(i) + all_arcs(k), 200);
    end
    arcHandles(k) = plot(r*cos(arc_t), r*sin(arc_t), 'Color', arcColors(k,:), 'LineWidth',2);
end

% === Spiral Trajectory (C^k for 0 <= k <= 9) ===
% Define fine-grained k values (continuous exponents)
k_vals = 0:0.05:9;
C      = A * exp(1i*2*pi*phi);
spiral = C.^k_vals;

% + spirial 
% % Normalize spiral the same way as dSet
spiral = spiral / sqrt(mean(abs(C.^marks).^2));

% Draw the spiral curve (background trajectory)
plot(real(spiral), imag(spiral), '-', ...
     'Color', 'k', 'LineWidth', 1.5, ...
     'DisplayName', 'C^k spiral');


% 5) Legend of arcs only
legend(arcHandles, arcLabels, 'Location','eastoutside', 'FontSize',9, 'NumColumns',1);

axis equal;
grid on;
box on;
xlabel('In-phase');
ylabel('Quadrature');
xlim([-5 5]);
ylim([-5 5]);
xticks([-5:1:5]);




%%
clear; clc;

% 1) Parameters
phi     = 1/13;                     % base phase
marks   = [0, 1, 3, 9];             % Golomb-ruler marks
A       = 0.95;
C       = A * exp(1i * 2 * pi * phi);
dSet    = C.^marks;                % 4-point constellation
dSet    = dSet / sqrt(mean(abs(dSet).^2));
M       = numel(dSet);

L        = 8;                      % symbols per frame
EbN0_dB  = 22;                     % demo SNR
numIter  = 2000;                   % Monte Carlo trials
tau_max  = 1e-6;                   % channel spread
delta_f  = 15e3;                   % frequency offset

% 2) Build equalized receive cloud (perfect CSIR)
eq_points = [];
for iter = 1:numIter
    idx     = randi([1 M], L, 1);
    d       = dSet(idx);
    delay_p = unifrnd(0, tau_max);
    IT      = randn + 1i * randn;
    h       = IT * exp(-2i * pi * delay_p * delta_f).^(0:L-1);
    rx      = awgn(h.' .* d, EbN0_dB, 'measured');
    eq      = rx ./ h.';      % perfect equalization
    eq_points = [eq_points; eq(:)];
end

% 3) Compute arcs: 6 CW + 6 CCW + full circle
angles = 2 * pi * phi * marks;
pair_list = nchoosek(1:M, 2);
Npairs = size(pair_list, 1);

d_cw  = zeros(Npairs,1);
d_ccw = zeros(Npairs,1);
k_cw  = zeros(Npairs,1);
k_ccw = zeros(Npairs,1);

for k = 1:Npairs
    i = pair_list(k,1);
    j = pair_list(k,2);
    % CW
    k_cw(k)  = mod(marks(j) - marks(i), (1/phi));
    d_cw(k)  = 2 * pi * k_cw(k) / (1/phi);
    % CCW
    k_ccw(k) = mod(marks(i) - marks(j), (1/phi));
    d_ccw(k) = 2 * pi * k_ccw(k) / (1/phi);
end
full_cw = 2*pi;
all_arcs = [d_cw; d_ccw; full_cw];
k_vals = [k_cw; k_ccw; (1/phi)];

arc_pairs = [pair_list; flip(pair_list,2); NaN NaN];
arcColors = hsv(numel(all_arcs));

% 4) 라벨링: "0→1 (2π·1/13)" 형식
arcLabels = cell(numel(all_arcs),1);
for k = 1:numel(all_arcs)
    if k <= Npairs || k <= 2*Npairs
        i = arc_pairs(k,1);
        j = arc_pairs(k,2);
        arcLabels{k} = sprintf('%d→%d (2\\pi·%d/13)', marks(i), marks(j), k_vals(k));
    else
        arcLabels{k} = 'full circle (2\pi)';
    end
end

% 5) Plot
figure;
hold on;

% a) Rx cloud
scatter(real(eq_points), imag(eq_points), 8, '.', 'MarkerEdgeAlpha',0.15, 'HandleVisibility','off');

% b) unit circle
theta = linspace(0,2*pi,360);
plot(cos(theta), sin(theta), 'k--','LineWidth',1,'HandleVisibility','off');

% c) Tx constellation
scatter(real(dSet), imag(dSet), 100, 'r+','LineWidth',1, 'HandleVisibility','off');
for i = 1:M
    plot([0 real(dSet(i))], [0 imag(dSet(i))], 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
    text(real(dSet(i))*1.15, imag(dSet(i))*1.15, sprintf('%d', marks(i)), ...
        'HorizontalAlignment','center', 'FontSize',12,'Color','g','FontWeight','bold');
end

% d) draw arcs with labels
arcHandles = gobjects(numel(all_arcs),1);
for k = 1:numel(all_arcs)
    r = 2 + 0.2 * k;
    i = arc_pairs(k,1);
    if isnan(i)
        arc_t = theta;
    else
        arc_t = linspace(angles(i), angles(i) + all_arcs(k), 200);
    end
    arcHandles(k) = plot(r * cos(arc_t), r * sin(arc_t), 'Color', arcColors(k,:), 'LineWidth',2);
end

% === Spiral Trajectory (C^k for 0 <= k <= 9) ===
% Define fine-grained k values (continuous exponents)
k_val = 0:0.01:9;
C      = A * exp(1i*2*pi*phi);
spiral = C.^k_val;

% + spirial 
% Normalize spiral the same way as dSet
spiral = spiral / sqrt(mean(abs(C.^marks).^2));

% Draw the spiral curve (background trajectory)
plot(real(spiral), imag(spiral), '-', ...
     'Color', 'k', 'LineWidth', 1.5, ...
     'DisplayName', 'C^k spiral');

% 6) Legend
legend(arcHandles, arcLabels, 'Location','eastoutside', 'FontSize',9, 'NumColumns',1);

axis equal;
grid on;
box on;
xlabel('In-phase');
ylabel('Quadrature');
xlim([-5 5]);
ylim([-5 5]);
xticks(-5:1:5);



%% Pairwise difference mapping
figure;
hold on;

colors = hsv(M * (M - 1) + 1);  % i ≠ j: M*(M-1)개 + i=j 하나
legend_handles = gobjects(M * (M - 1) + 1, 1);
label_idx = 1;
v = 13;

for i = 1:M
    for j = 1:M
        
        diff_pt = dSet(i) / dSet(j);  % C^{m_i - m_j}
        diff_k = mod(marks(j) - marks(i), v);  % 0~12

        if i == j
            if i == 1  
                legend_handles(label_idx) = scatter(real(diff_pt), imag(diff_pt), ...
                    80, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3], ...
                    'DisplayName', '(2\pi·0/13)');
                label_idx = label_idx + 1;
            end
        else
            % 일반적인 pairwise 차이
            legend_handles(label_idx) = scatter(real(diff_pt), imag(diff_pt), ...
                60, 'filled', 'MarkerFaceColor', colors(label_idx,:), ...
                'DisplayName', sprintf('%d→%d (2\\pi·%d/13)', marks(i), marks(j), diff_k));
            label_idx = label_idx + 1;
        end
    end
end

% 배경 unit circle
theta = linspace(0, 2*pi, 360);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1, 'HandleVisibility','off');
plot(0, 0, 'k+', 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility','off');

% 설정
axis equal;
grid on;
box on;
xlabel('In-phase');
ylabel('Quadrature');
xlim([-2 2]);
ylim([-2 2]);
xticks(-5:1:5);
yticks(-5:1:5);

% Legend
legend(legend_handles, 'Location','eastoutside', 'FontSize', 9);
