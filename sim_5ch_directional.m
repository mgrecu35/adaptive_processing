function [X, meta] = sim_5ch_directional(range_m, sim_prms, clutter_state)
%SIM_5CH_DIRECTIONAL  5-channel pulse-train simulation with direction-dependent
%                     beam patterns and Rice-compound distributed clutter.
%
% Clutter amplitude uses a compound Rice model:
%
%   X_clut(t) = sqrt( m^2/(m^2+1) * X1^2  +  1/(m^2+1) * X2(t)^2 )
%
%   X1  ~ Weibull(k_w, lambda_w)  — persistent component, fixed per scatterer
%   X2  ~ N(0, sigma^2)           — fluctuating component, redrawn each scan
%   sigma = lambda_w * sqrt( Gamma(1 + 2/k_w) )  — preserves total power
%   m   — Rice parameter: large m → persistent, m→0 → fully fluctuating
%
% The persistent component X1 and the scatterer directions are held in
% CLUTTER_STATE across successive calls (scans).  Pass [] on the first call
% and the returned meta.clutter_state on subsequent calls.
%
% INPUTS
%   range_m       - range of the cell under test [m]
%   sim_prms      - parameter struct from build_sim_prms_dir()
%   clutter_state - [] for first scan, or meta.clutter_state from prior scan
%
% OUTPUTS
%   X         - [N_ch x N_pulses] complex IQ matrix (one row per channel)
%   meta      - struct with individual contributions and geometry
%               meta.clutter_state  — pass to next call to maintain persistence

arguments
    range_m       (1,1) double = 5000
    sim_prms            struct = build_sim_prms_dir()
    clutter_state              = []     % [] → first scan; struct → subsequent scans
end

%% ── 0. Unpack parameters ────────────────────────────────────────────────────
lambda      = sim_prms.lambda;
N_ch        = sim_prms.N_ch;          % 5
N_pulses    = sim_prms.N_pulses;      % 64 | 128 | 256
PRI         = sim_prms.PRI;
v_target    = sim_prms.v_target;      % radial velocity [m/s]
SNR_dB      = sim_prms.SNR_dB;
CNR_dB      = sim_prms.CNR_dB;
noise_pwr   = sim_prms.noise_pwr;
n_clutter   = sim_prms.n_clutter_patches;
weibull_k   = sim_prms.weibull_k;
rice_m      = sim_prms.rice_m;        % Rice parameter (persistent/fluctuating ratio)
clutter_az  = sim_prms.clutter_az;   % nominal clutter az [rad]
clutter_el  = sim_prms.clutter_el;   % nominal clutter el [rad]
clutter_az_spread = sim_prms.clutter_az_spread;  % angular std [rad]
clutter_el_spread = sim_prms.clutter_el_spread;
tgt_az      = sim_prms.target_az;
tgt_el      = sim_prms.target_el;
patch_pos   = sim_prms.patch_pos;    % [3 x 41]
sub_arrays  = sim_prms.sub_arrays;   % [5 x 41] logical

%% ── 1. Subarray geometry ────────────────────────────────────────────────────
% Phase centres (unweighted mean of patch positions per subarray)
rx_pos = zeros(3, N_ch);
for m = 1:N_ch
    rx_pos(:, m) = mean(patch_pos(:, sub_arrays(m,:)), 2);
end

% Effective aperture diameter per subarray (RMS spread of patch positions)
D = zeros(1, N_ch);
for m = 1:N_ch
    dp = patch_pos(:, sub_arrays(m,:)) - rx_pos(:, m);
    D(m) = 2 * sqrt(mean(sum(dp.^2, 1)));
end
% Channel 5 has D=0 (single patch) → treat as omnidirectional below

% Gaussian beamwidth parameter σ [rad]: σ ≈ λ/(π·D)
% For D=0 set σ→∞ (gain = 1 for all angles)
sigma_beam = zeros(1, N_ch);
for m = 1:N_ch
    if D(m) > 0
        sigma_beam(m) = lambda / (pi * D(m));
    else
        sigma_beam(m) = Inf;   % omnidirectional
    end
end

%% ── 2. Helper: unit vector, steering vector, beam gain ─────────────────────
% Unit vector for (az, el)
uvec = @(az, el) [cos(el).*cos(az); cos(el).*sin(az); sin(el)];

% Steering vector for an array of phase centres rx_pos [3 x M]
steer = @(az, el) exp(1j * (2*pi/lambda) * (rx_pos.' * uvec(az, el)));

% Gaussian beam gain for channel m at angular offset δθ from boresight
beam_gain = @(delta_theta, m) exp(-delta_theta.^2 / (2 * sigma_beam(m)^2));

% Angular offset between two (az,el) directions [scalar]
ang_offset = @(az1, el1, az2, el2) acos( ...
    max(-1, min(1, dot(uvec(az1,el1), uvec(az2,el2)) )));

%% ── 3. Noise power and signal scaling ──────────────────────────────────────
% Reference gain: coherent sum of N_patches in main subarray
g_main = sqrt(10);   % subarrays 1-4
% (aux channel has gain 1 at boresight by the Gaussian model too)

SNR_lin = 10^(SNR_dB/10);
CNR_lin = 10^(CNR_dB/10);

% Target amplitude: set so that main-channel output SNR = SNR_dB
% |alpha|^2 * g_main^2 = SNR_lin * noise_pwr
alpha_mag = sqrt(SNR_lin * noise_pwr) / g_main;

% Clutter scale: total clutter power in main channel at boresight = CNR_lin * noise_pwr
% Distributed over n_clutter patches.
%
% Second moment of Weibull(k_w, lambda_w):  E[X1^2] = lambda_w^2 * Gamma(1+2/k_w)
% sigma^2 is set equal to E[X1^2] to preserve total power across the compound model.
% E[X^2] = m^2/(m^2+1)*E[X1^2] + 1/(m^2+1)*sigma^2 = E[X1^2]  ✓ (independent of m)
clutter_pwr_total = CNR_lin * noise_pwr;
weibull_mean2 = gamma(1 + 2/weibull_k);        % E[a^2] for unit-scale Weibull
lambda_w      = sqrt(clutter_pwr_total / (n_clutter * weibull_mean2));
sigma_fluct   = lambda_w * sqrt(weibull_mean2); % = lambda_w*sqrt(Gamma(1+2/k_w))

%% ── 4. Clutter scatterer geometry and persistence ───────────────────────────
%
% On the FIRST scan (clutter_state is empty):
%   - Fix scatterer directions and Weibull (persistent) amplitudes X1.
%   - Draw initial random phases.
%
% On SUBSEQUENT scans:
%   - Reuse directions and X1 from clutter_state.
%   - Redraw the fluctuating component X2 ~ N(0, sigma_fluct^2).
%   - Combine via the Rice compound model.

if isempty(clutter_state)
    % ── First scan: initialise persistent state ──────────────────────────
    rng(sim_prms.rng_seed);

    % Scatterer directions — fixed for all scans
    az_c = clutter_az + clutter_az_spread * randn(1, n_clutter);
    el_c = clutter_el + clutter_el_spread * randn(1, n_clutter);

    % Weibull persistent amplitudes X1 — fixed for all scans
    u    = rand(1, n_clutter);
    X1_c = lambda_w * (-log(1 - u)).^(1/weibull_k);   % [1 x n_clutter]

    % Persistent phase — fixed for all scans.
    % This is the key quantity that encodes inter-channel spatial coherence:
    % x_clut(m) = sum_i amp_i * exp(j*phase_i) * G(m,i) * exp(j*2pi/lambda * u_i'*p_m)
    % The spatial covariance E[x_clut * x_clut^H] is nonzero off-diagonal
    % ONLY because phase_c is consistent across scans.  Redrawing phase_c
    % every scan causes the cross-terms to average to zero → diagonal R.
    phase_c = 2*pi * rand(1, n_clutter);

else
    % ── Subsequent scans: restore persistent state ────────────────────────
    az_c    = clutter_state.az_c;
    el_c    = clutter_state.el_c;
    X1_c    = clutter_state.X1_c;
    phase_c = clutter_state.phase_c;   % ← restored, not redrawn
end

% Fluctuating component X2 ~ N(0, sigma_fluct^2) — redrawn every scan.
% This modulates the amplitude scan-to-scan but does NOT scramble the phase,
% so the spatial correlation structure of x_clut is preserved.
X2_c  = sigma_fluct * randn(1, n_clutter);    % [1 x n_clutter]

% Rice compound amplitude:
%   amp_c = sqrt( m^2/(m^2+1)*X1^2 + 1/(m^2+1)*X2^2 )
m2    = rice_m^2;
amp_c = sqrt( (m2 * X1_c.^2 + X2_c.^2) / (m2 + 1) );  % [1 x n_clutter]

%% ── 5. Per-channel beam gain for every scatterer ───────────────────────────
% G_clutter(m, i) = Gaussian gain of channel m toward scatterer i
G_clut = zeros(N_ch, n_clutter);
for m = 1:N_ch
    for i = 1:n_clutter
        dtheta = ang_offset(az_c(i), el_c(i), tgt_az, tgt_el);
        % boresight is the array pointing direction (target direction here)
        G_clut(m, i) = beam_gain(dtheta, m);
    end
end

% Steering vectors for all clutter scatterers [N_ch x n_clutter]
S_clut = zeros(N_ch, n_clutter);
for i = 1:n_clutter
    S_clut(:, i) = steer(az_c(i), el_c(i));
end

% Clutter spatial snapshot (no Doppler — stationary ground clutter)
% x_clut(m) = sum_i  amp_i * exp(j*phase_i) * G(m,i) * s_m(az_i, el_i)
clut_weights = amp_c .* exp(1j * phase_c);            % [1 x n_clutter]
x_clut = (G_clut .* S_clut) * clut_weights.';         % [N_ch x 1]

%% ── 6. Target steering vector and amplitude ─────────────────────────────────
s_tgt   = steer(tgt_az, tgt_el);                      % [N_ch x 1]

% Per-channel target gain: beam gain of each subarray toward the target
% Target is ON boresight → G_tgt(m) = 1 for all m (by definition)
% (If target were off-boresight you would compute ang_offset here)
g_tgt   = ones(N_ch, 1);
g_tgt(5) = beam_gain(ang_offset(tgt_az, tgt_el, tgt_az, tgt_el), 5); % = 1

alpha   = alpha_mag * exp(1j * 2*pi*rand());

%% ── 7. Build pulse train ────────────────────────────────────────────────────
f_d   = 2 * v_target / lambda;          % Doppler frequency [Hz]
n_vec = (0 : N_pulses-1).';             % pulse index column vector

% Doppler phase per pulse [N_pulses x 1]
doppler_phase = exp(1j * 2*pi * f_d * n_vec * PRI);

% Target: [N_ch x N_pulses]  — same spatial snapshot, rotating phase
X_tgt  = (alpha * (g_tgt .* s_tgt)) * doppler_phase.';

% Clutter: [N_ch x N_pulses]  — stationary (same snapshot every pulse)
%   Add slow internal clutter motion as small pulse-to-pulse phase jitter
if sim_prms.clutter_doppler_std > 0
    clut_phase_jitter = exp(1j * sim_prms.clutter_doppler_std * randn(1, N_pulses));
    X_clut = x_clut * clut_phase_jitter;
else
    X_clut = repmat(x_clut, [1, N_pulses]);
end

% Noise: [N_ch x N_pulses]
X_noise = sqrt(noise_pwr/2) * (randn(N_ch, N_pulses) + 1j*randn(N_ch, N_pulses));

% Total received matrix
X = X_tgt + X_clut + X_noise;

%% ── 8. Package metadata ─────────────────────────────────────────────────────
meta.rx_pos        = rx_pos;
meta.D             = D;
meta.sigma_beam    = sigma_beam;
meta.G_clut        = G_clut;      % [N_ch x n_clutter] beam gains
meta.az_c          = az_c;
meta.el_c          = el_c;
meta.amp_c         = amp_c;       % compound Rice amplitudes this scan
meta.X1_c          = X1_c;       % persistent Weibull amplitudes
meta.X2_c          = X2_c;       % fluctuating component this scan
meta.x_clut        = x_clut;
meta.X_tgt         = X_tgt;
meta.X_clut        = X_clut;
meta.X_noise       = X_noise;
meta.f_d           = f_d;
meta.s_tgt         = s_tgt;
meta.lambda_w      = lambda_w;
meta.sigma_fluct   = sigma_fluct;
meta.rice_m        = rice_m;

% Persistent state — pass to next call to maintain scan-to-scan continuity
meta.clutter_state.az_c    = az_c;
meta.clutter_state.el_c    = el_c;
meta.clutter_state.X1_c    = X1_c;
meta.clutter_state.phase_c = phase_c;   % ← essential for spatial coherence

end % sim_5ch_directional

function [X, meta] = sim_5ch_directional_old(range_m, sim_prms, clutter_state)
%SIM_5CH_DIRECTIONAL  5-channel pulse-train simulation with direction-dependent
%                     beam patterns and Rice-compound distributed clutter.
%
% Clutter amplitude uses a compound Rice model:
%
%   X_clut(t) = sqrt( m^2/(m^2+1) * X1^2  +  1/(m^2+1) * X2(t)^2 )
%
%   X1  ~ Weibull(k_w, lambda_w)  — persistent component, fixed per scatterer
%   X2  ~ N(0, sigma^2)           — fluctuating component, redrawn each scan
%   sigma = lambda_w * sqrt( Gamma(1 + 2/k_w) )  — preserves total power
%   m   — Rice parameter: large m → persistent, m→0 → fully fluctuating
%
% The persistent component X1 and the scatterer directions are held in
% CLUTTER_STATE across successive calls (scans).  Pass [] on the first call
% and the returned meta.clutter_state on subsequent calls.
%
% INPUTS
%   range_m       - range of the cell under test [m]
%   sim_prms      - parameter struct from build_sim_prms_dir()
%   clutter_state - [] for first scan, or meta.clutter_state from prior scan
%
% OUTPUTS
%   X         - [N_ch x N_pulses] complex IQ matrix (one row per channel)
%   meta      - struct with individual contributions and geometry
%               meta.clutter_state  — pass to next call to maintain persistence

arguments
    range_m       (1,1) double = 5000
    sim_prms            struct = build_sim_prms_dir()
    clutter_state              = []     % [] → first scan; struct → subsequent scans
end

%% ── 0. Unpack parameters ────────────────────────────────────────────────────
lambda      = sim_prms.lambda;
N_ch        = sim_prms.N_ch;          % 5
N_pulses    = sim_prms.N_pulses;      % 64 | 128 | 256
PRI         = sim_prms.PRI;
v_target    = sim_prms.v_target;      % radial velocity [m/s]
SNR_dB      = sim_prms.SNR_dB;
CNR_dB      = sim_prms.CNR_dB;
noise_pwr   = sim_prms.noise_pwr;
n_clutter   = sim_prms.n_clutter_patches;
weibull_k   = sim_prms.weibull_k;
rice_m      = sim_prms.rice_m;        % Rice parameter (persistent/fluctuating ratio)
clutter_az  = sim_prms.clutter_az;   % nominal clutter az [rad]
clutter_el  = sim_prms.clutter_el;   % nominal clutter el [rad]
clutter_az_spread = sim_prms.clutter_az_spread;  % angular std [rad]
clutter_el_spread = sim_prms.clutter_el_spread;
tgt_az      = sim_prms.target_az;
tgt_el      = sim_prms.target_el;
patch_pos   = sim_prms.patch_pos;    % [3 x 41]
sub_arrays  = sim_prms.sub_arrays;   % [5 x 41] logical

%% ── 1. Subarray geometry ────────────────────────────────────────────────────
% Phase centres (unweighted mean of patch positions per subarray)
rx_pos = zeros(3, N_ch);
for m = 1:N_ch
    rx_pos(:, m) = mean(patch_pos(:, sub_arrays(m,:)), 2);
end

% Effective aperture diameter per subarray (RMS spread of patch positions)
D = zeros(1, N_ch);
for m = 1:N_ch
    dp = patch_pos(:, sub_arrays(m,:)) - rx_pos(:, m);
    D(m) = 2 * sqrt(mean(sum(dp.^2, 1)));
end
% Channel 5 has D=0 (single patch) → treat as omnidirectional below

% Gaussian beamwidth parameter σ [rad]: σ ≈ λ/(π·D)
% For D=0 set σ→∞ (gain = 1 for all angles)
sigma_beam = zeros(1, N_ch);
for m = 1:N_ch
    if D(m) > 0
        sigma_beam(m) = lambda / (pi * D(m));
    else
        sigma_beam(m) = Inf;   % omnidirectional
    end
end

%% ── 2. Helper: unit vector, steering vector, beam gain ─────────────────────
% Unit vector for (az, el)
uvec = @(az, el) [cos(el).*cos(az); cos(el).*sin(az); sin(el)];

% Steering vector for an array of phase centres rx_pos [3 x M]
steer = @(az, el) exp(1j * (2*pi/lambda) * (rx_pos.' * uvec(az, el)));

% Gaussian beam gain for channel m at angular offset δθ from boresight
beam_gain = @(delta_theta, m) exp(-delta_theta.^2 / (2 * sigma_beam(m)^2));

% Angular offset between two (az,el) directions [scalar]
ang_offset = @(az1, el1, az2, el2) acos( ...
    max(-1, min(1, dot(uvec(az1,el1), uvec(az2,el2)) )));

%% ── 3. Noise power and signal scaling ──────────────────────────────────────
% Reference gain: coherent sum of N_patches in main subarray
g_main = sqrt(10);   % subarrays 1-4
% (aux channel has gain 1 at boresight by the Gaussian model too)

SNR_lin = 10^(SNR_dB/10);
CNR_lin = 10^(CNR_dB/10);

% Target amplitude: set so that main-channel output SNR = SNR_dB
% |alpha|^2 * g_main^2 = SNR_lin * noise_pwr
alpha_mag = sqrt(SNR_lin * noise_pwr) / g_main;

% Clutter scale: total clutter power in main channel at boresight = CNR_lin * noise_pwr
% Distributed over n_clutter patches.
%
% Second moment of Weibull(k_w, lambda_w):  E[X1^2] = lambda_w^2 * Gamma(1+2/k_w)
% sigma^2 is set equal to E[X1^2] to preserve total power across the compound model.
% E[X^2] = m^2/(m^2+1)*E[X1^2] + 1/(m^2+1)*sigma^2 = E[X1^2]  ✓ (independent of m)
clutter_pwr_total = CNR_lin * noise_pwr;
weibull_mean2 = gamma(1 + 2/weibull_k);        % E[a^2] for unit-scale Weibull
lambda_w      = sqrt(clutter_pwr_total / (n_clutter * weibull_mean2));
sigma_fluct   = lambda_w * sqrt(weibull_mean2); % = lambda_w*sqrt(Gamma(1+2/k_w))

%% ── 4. Clutter scatterer geometry and persistence ───────────────────────────
%
% On the FIRST scan (clutter_state is empty):
%   - Fix scatterer directions and Weibull (persistent) amplitudes X1.
%   - Draw initial random phases.
%
% On SUBSEQUENT scans:
%   - Reuse directions and X1 from clutter_state.
%   - Redraw the fluctuating component X2 ~ N(0, sigma_fluct^2).
%   - Combine via the Rice compound model.

if isempty(clutter_state)
    % ── First scan: initialise persistent state ──────────────────────────
    rng(sim_prms.rng_seed);

    % Scatterer directions — fixed for all scans
    az_c = clutter_az + clutter_az_spread * randn(1, n_clutter);
    el_c = clutter_el + clutter_el_spread * randn(1, n_clutter);

    % Weibull persistent amplitudes X1 — fixed for all scans
    u    = rand(1, n_clutter);
    X1_c = lambda_w * (-log(1 - u)).^(1/weibull_k);   % [1 x n_clutter]

else
    % ── Subsequent scans: restore persistent state ────────────────────────
    az_c = clutter_state.az_c;
    el_c = clutter_state.el_c;
    X1_c = clutter_state.X1_c;
end

% Fluctuating component X2 ~ N(0, sigma_fluct^2) — redrawn every scan
X2_c = sigma_fluct * randn(1, n_clutter);             % [1 x n_clutter]

% Rice compound amplitude:
%   amp_c = sqrt( m^2/(m^2+1)*X1^2 + 1/(m^2+1)*X2^2 )
m2   = rice_m^2;
amp_c = sqrt( (m2 * X1_c.^2 + X2_c.^2) / (m2 + 1) );  % [1 x n_clutter]

% Random phase — redrawn every scan (X2 carries both amplitude and phase
% fluctuation; the overall phase is still uniformly distributed)
phase_c = 2*pi * rand(1, n_clutter);

%% ── 5. Per-channel beam gain for every scatterer ───────────────────────────
% G_clutter(m, i) = Gaussian gain of channel m toward scatterer i
G_clut = zeros(N_ch, n_clutter);
for m = 1:N_ch
    for i = 1:n_clutter
        dtheta = ang_offset(az_c(i), el_c(i), tgt_az, tgt_el);
        % boresight is the array pointing direction (target direction here)
        G_clut(m, i) = beam_gain(dtheta, m);
    end
end

% Steering vectors for all clutter scatterers [N_ch x n_clutter]
S_clut = zeros(N_ch, n_clutter);
for i = 1:n_clutter
    S_clut(:, i) = steer(az_c(i), el_c(i));
end

% Clutter spatial snapshot (no Doppler — stationary ground clutter)
% x_clut(m) = sum_i  amp_i * exp(j*phase_i) * G(m,i) * s_m(az_i, el_i)
clut_weights = amp_c .* exp(1j * phase_c);            % [1 x n_clutter]
x_clut = (G_clut .* S_clut) * clut_weights.';         % [N_ch x 1]

%% ── 6. Target steering vector and amplitude ─────────────────────────────────
s_tgt   = steer(tgt_az, tgt_el);                      % [N_ch x 1]

% Per-channel target gain: beam gain of each subarray toward the target
% Target is ON boresight → G_tgt(m) = 1 for all m (by definition)
% (If target were off-boresight you would compute ang_offset here)
g_tgt   = ones(N_ch, 1);
g_tgt(5) = beam_gain(ang_offset(tgt_az, tgt_el, tgt_az, tgt_el), 5); % = 1

alpha   = alpha_mag * exp(1j * 2*pi*rand());

%% ── 7. Build pulse train ────────────────────────────────────────────────────
f_d   = 2 * v_target / lambda;          % Doppler frequency [Hz]
n_vec = (0 : N_pulses-1).';             % pulse index column vector

% Doppler phase per pulse [N_pulses x 1]
doppler_phase = exp(1j * 2*pi * f_d * n_vec * PRI);

% Target: [N_ch x N_pulses]  — same spatial snapshot, rotating phase
X_tgt  = (alpha * (g_tgt .* s_tgt)) * doppler_phase.';

% Clutter: [N_ch x N_pulses]  — stationary (same snapshot every pulse)
%   Add slow internal clutter motion as small pulse-to-pulse phase jitter
if sim_prms.clutter_doppler_std > 0
    clut_phase_jitter = exp(1j * sim_prms.clutter_doppler_std * randn(1, N_pulses));
    X_clut = x_clut * clut_phase_jitter;
else
    X_clut = repmat(x_clut, [1, N_pulses]);
end

% Noise: [N_ch x N_pulses]
X_noise = sqrt(noise_pwr/2) * (randn(N_ch, N_pulses) + 1j*randn(N_ch, N_pulses));

% Total received matrix
X = X_tgt + X_clut + X_noise;

%% ── 8. Package metadata ─────────────────────────────────────────────────────
meta.rx_pos        = rx_pos;
meta.D             = D;
meta.sigma_beam    = sigma_beam;
meta.G_clut        = G_clut;      % [N_ch x n_clutter] beam gains
meta.az_c          = az_c;
meta.el_c          = el_c;
meta.amp_c         = amp_c;       % compound Rice amplitudes this scan
meta.X1_c          = X1_c;       % persistent Weibull amplitudes
meta.X2_c          = X2_c;       % fluctuating component this scan
meta.x_clut        = x_clut;
meta.X_tgt         = X_tgt;
meta.X_clut        = X_clut;
meta.X_noise       = X_noise;
meta.f_d           = f_d;
meta.s_tgt         = s_tgt;
meta.lambda_w      = lambda_w;
meta.sigma_fluct   = sigma_fluct;
meta.rice_m        = rice_m;

% Persistent state — pass to next call to maintain scan-to-scan continuity
meta.clutter_state.az_c = az_c;
meta.clutter_state.el_c = el_c;
meta.clutter_state.X1_c = X1_c;

end % sim_5ch_directional

