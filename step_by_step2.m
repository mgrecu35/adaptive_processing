prms            = build_sim_prms_dir();
prms.N_pulses   = 64;
prms.v_target   = 6;       % keep target away from zero Doppler
prms.SNR_dB     = 15;
prms.CNR_dB     = 15;
N_scans         = 50;
prms.rice_m     = 3;
N               = prms.N_pulses;

% --- accumulate covariance from zero-Doppler bin, channels 1:4 only ---
R_hat   = zeros(4, 4);
cs      = [];
X_dop_last = [];

for t = 1:N_scans
    [X, meta_t]  = sim_5ch_directional(5000, prms, cs);
    cs           = meta_t.clutter_state;
    X_dop        = calc_doppler(X, prms);   % [5 x N]

    % zero-Doppler bin (bin 1, pre-fftshift), channels 1:4
    x_clut       = X_dop(1:4, 1);          % [4 x 1]
    R_hat        = R_hat + (x_clut * x_clut') / N_scans;

    X_dop_last   = X_dop;                  % keep last scan for plotting
end

% --- diagonal loading: 10 dB below mean eigenvalue ---
delta   = trace(R_hat) / 4 * 10^(-10/10);
R_reg   = R_hat + delta * eye(4);

% --- MVDR weight vector ---
s       = ones(4, 1);                      % boresight steering, channels 1:4
R_inv_s = R_reg \ s;
w       = R_inv_s / (s' * R_inv_s);       % [4 x 1], distortionless constraint

% --- apply filter across ALL Doppler bins of last scan ---
X_mvdr  = w' * X_dop_last(1:4, :);        % [1 x N]

% --- conventional sum beam for comparison ---
X_sum   = mean(X_dop_last(1:4, :), 1);    % [1 x N]

% --- plot ---
f_axis  = (-N/2 : N/2-1) / (N * prms.PRI);
v_axis  = f_axis * prms.lambda / 2;

X_mvdr_shift = fftshift(X_mvdr);
X_sum_shift  = fftshift(X_sum);
X_dop_shift  = fftshift(X_dop_last, 2);

figure('Name','MVDR vs Sum Beam','NumberTitle','off');
hold on;
colors = lines(4);
for m = 1:4
    plot(v_axis, 20*log10(abs(X_dop_shift(m,:))+eps), ...
        'Color', colors(m,:), 'LineWidth', 0.8);
end
plot(v_axis, 20*log10(abs(X_sum_shift)+eps),  'k--', 'LineWidth', 1.5);
plot(v_axis, 20*log10(abs(X_mvdr_shift)+eps), 'r-',  'LineWidth', 2.0);
xlabel('Velocity [m/s]'); ylabel('Power [dB]'); grid on;
legend('Ch1','Ch2','Ch3','Ch4','Sum beam','MVDR');
title('MVDR adaptive filter vs conventional sum beam');
