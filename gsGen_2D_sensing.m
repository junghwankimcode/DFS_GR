function [GS, Y_ideal, Y_noisy, y_target, channel, R_true, V_true] = gsGen_2D_sensing(SNR, S, P1, P2)
    % TF-domain sensing channel generator (your structure kept)
    carrier_freq = 28e9;
    c_speed      = physconst("lightspeed");
    delta_f      = 480e3;
    T            = 1 / delta_f;

    R_limit = 200;     % meters
    V_limit = 20;      % m/s

    % Rician K-factor
    K_dB  = 6;
    K_lin = 10^(K_dB/10);

    % ranges (LoS is the closest)
    ranges = sort(unifrnd(0, R_limit, [1, S]));
    R_true = ranges(1);

    channel = zeros(P1, P2);
    n = (0:P1-1).';
    m = (0:P2-1);

    GS = cell(1, S);
    V_true = 0;

    for ii = 1:S
        if ii == 1
            V_true   = unifrnd(-V_limit, V_limit);
            target_v = V_true;

            if S == 1
                gain_val = exp(1i*2*pi*rand(1));
            else
                gain_val = sqrt(K_lin/(K_lin+1)) * exp(1i*2*pi*rand(1));
            end
        else
            % carrier_freq * v_max * cos(2 * pi * rand(1)) / c_speed;
            target_v = unifrnd(-V_limit, V_limit);
            raw_gain = (randn(1) + 1i*randn(1)) / sqrt(2);
            gain_val = sqrt(1/(K_lin+1)) * raw_gain / sqrt(S-1);
        end

        % round-trip delay & Doppler
        delay_p = 2 * ranges(ii) / c_speed;
        angle2  = delay_p * delta_f;

        f_D     = (2 * (target_v) * carrier_freq) / c_speed;
        angle1  = f_D * T;

        GS{ii}.gain   = gain_val;
        GS{ii}.angle1 = angle1;
        GS{ii}.angle2 = angle2;

        GS{ii}.Form = gain_val ...
            * exp(1i * 2*pi * angle1 * n) ...
            * exp(-1i * 2*pi * angle2 * m);

        channel = channel + GS{ii}.Form;
    end

    % output
    X       = ones(P1, P2);
    Y_ideal = X .* channel;
    Y_noisy = awgn(Y_ideal, SNR, 'measured');

    % keep your legacy downsampling outputs (not used in main loop)
    y_target = Y_noisy(2:2:end, 2:2:end);
    Y_ideal  = Y_ideal(2:2:end, 2:2:end);
end
% function [GS, Y_ideal, Y_noisy, y_target, channel, R_true, V_true, info] = gsGen_2D_sensing(SNR, S, P1, P2)
% 
%     carrier_freq = 28e9;
%     c_speed      = physconst("lightspeed");
%     delta_f      = 480e3;
%     T            = 1 / delta_f;
% 
%     R_limit = 200;
%     V_limit = 20;
% 
%     % Rician K-factor
%     K_dB  = 6;
%     K_lin = 10^(K_dB/10);
% 
%     % path physical params
%     ranges = sort(unifrnd(0, R_limit, [1, S]));  % LoS is shortest
%     v_all  = unifrnd(-V_limit, V_limit, [1, S]);
% 
%     channel = zeros(P1, P2);
%     n = (0:P1-1).';
%     m = (0:P2-1);
% 
%     GS = cell(1, S);
% 
%     alpha = zeros(1,S);
% 
%     for ii = 1:S
%         if ii == 1
%             % LoS (dominant component)
%             if S == 1
%                 alpha(ii) = exp(1i*2*pi*rand(1));
%             else
%                 alpha(ii) = sqrt(K_lin/(K_lin+1)) * exp(1i*2*pi*rand(1));
%             end
%             gain_val = alpha(ii);
%         else
%             % NLoS (scattering)
%             raw_gain = (randn(1) + 1i*randn(1)) / sqrt(2);
%             alpha(ii) = sqrt(1/(K_lin+1)) * raw_gain / sqrt(S-1);
%             gain_val  = alpha(ii);
%         end
% 
%         delay_p = 2 * ranges(ii) / c_speed;   % round-trip
%         angle2  = delay_p * delta_f;
% 
%         f_D     = (2 * v_all(ii) * carrier_freq) / c_speed;
%         angle1  = f_D * T;
% 
%         GS{ii}.gain   = gain_val;
%         GS{ii}.angle1 = angle1;
%         GS{ii}.angle2 = angle2;
% 
%         GS{ii}.Form = gain_val ...
%             * exp(1i * 2*pi * angle1 * n) ...
%             * exp(-1i * 2*pi * angle2 * m);
% 
%         channel = channel + GS{ii}.Form;
%     end
% 
%     % === IMPORTANT: GT = strongest path (matches "global peak" estimator) ===
%     [~, idx_str] = max(abs(alpha));
%     R_true = ranges(idx_str);
%     V_true = v_all(idx_str);
% 
%     info.ranges = ranges;
%     info.v_all  = v_all;
%     info.alpha  = alpha;
%     info.idx_str = idx_str;
% 
%     % output
%     X       = ones(P1, P2);
%     Y_ideal = X .* channel;
%     Y_noisy = awgn(Y_ideal, SNR, 'measured');
% 
%     y_target = Y_noisy(2:2:end, 2:2:end);
%     Y_ideal  = Y_ideal(2:2:end, 2:2:end);
% end
