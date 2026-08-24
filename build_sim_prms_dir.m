%% ═══════════════════════════════════════════════════════════════════════════
function prms = build_sim_prms_dir()
%BUILD_SIM_PRMS_DIR  Default simulation parameters.

%% Radar / waveform
prms.lambda     = 3e8 / 3.4e9;   % wavelength [m]  (~0.0882 m)
prms.N_ch       = 5;
prms.N_pulses   = 64;             % {64, 128, 256}
prms.PRI        = 1e-3;           % [s]  → PRF = 1 kHz

%% Target
prms.target_az  = 0;             % [rad] on boresight
prms.target_el  = 0;             % [rad]
prms.v_target   = 30;            % radial velocity [m/s]
prms.SNR_dB     = 15;            % desired per-element SNR [dB]

%% Clutter
prms.clutter_az         = 0;     % nominal clutter az (same as beam centre)
prms.clutter_el         = 0;
prms.clutter_az_spread  = deg2rad(30);  % wide spread — distributed ground clutter
prms.clutter_el_spread  = deg2rad(5);
prms.n_clutter_patches  = 200;          % number of independent scatterers
prms.CNR_dB             = 30;           % desired total CNR in main channel [dB]
prms.weibull_k          = 0.8;          % shape: <1 → heavier tail than Rayleigh
prms.rice_m             = 3.0;          % Rice parameter: large → persistent clutter
%                                       % m=0 → Rayleigh (fully fluctuating)
%                                       % m>>1 → near-deterministic (frozen)
%                                       % scan-to-scan amplitude correlation ≈ m^2/(m^2+1)
prms.clutter_doppler_std = 0.05;        % pulse-to-pulse phase jitter std [rad]

%% Noise
prms.noise_pwr  = 1.0;          % normalised thermal noise power per channel

%% Array geometry — extracted from set_rdr_cfgs
%  patch_pos: [3 x 41], columns are patch positions in metres
%  Coordinates: row 1 = x (boresight), row 2 = y (lateral), row 3 = z (vertical)
patch_pos_mm = [ ...
       0,         0,          0; ...  1
       0,   054.870,   -022.728; ...  2
       0,   054.870,    022.728; ...  3
       0,   022.728,    054.870; ...  4
       0,  -022.728,    054.870; ...  5
       0,  -054.870,    022.728; ...  6
       0,  -054.870,   -022.728; ...  7
       0,  -022.728,   -054.870; ...  8
       0,   022.728,   -054.870; ...  9
       0,   094.236,          0; ... 10
       0,   066.635,    066.635; ... 11
       0,         0,    094.236; ... 12
       0,  -066.635,    066.635; ... 13
       0,  -094.236,          0; ... 14
       0,  -066.635,   -066.635; ... 15
       0,         0,   -094.236; ... 16
       0,   066.635,   -066.635; ... 17
       0,   107.624,   -044.579; ... 18
       0,   107.624,    044.579; ... 19
       0,   044.579,    107.624; ... 20
       0,  -044.579,    107.624; ... 21
       0,  -107.624,    044.579; ... 22
       0,  -107.624,   -044.579; ... 23
       0,  -044.579,   -107.624; ... 24
       0,   044.579,   -107.624; ... 25
       0,   139.691,          0; ... 26
       0,   098.776,    098.776; ... 27
       0,         0,    139.691; ... 28
       0,  -098.776,    098.776; ... 29
       0,  -139.691,          0; ... 30
       0,  -098.776,   -098.776; ... 31
       0,         0,   -139.691; ... 32
       0,   098.776,   -098.776; ... 33
       0,   149.618,   -061.974; ... 34
       0,   149.618,    061.974; ... 35
       0,   061.974,    149.618; ... 36
       0,  -061.974,    149.618; ... 37
       0,  -149.618,    061.974; ... 38
       0,  -149.618,   -061.974; ... 39
       0,  -061.974,   -149.618; ... 40
       0,   061.974,   -149.618; ... 41
       ].' * 1e-3;

prms.patch_pos = patch_pos_mm;

% Subarray membership — logical [5 x 41]
sub_arrays = false(5, 41);
sub_arrays(1, [3,4,11,12,19,20,26,27,35,36]) = true;
sub_arrays(2, [5,6,13,14,21,22,28,29,37,38]) = true;
sub_arrays(3, [7,8,15,16,23,24,30,31,39,40]) = true;
sub_arrays(4, [2,9,10,17,18,25,32,33,34,41]) = true;
sub_arrays(5, [1])                            = true;
prms.sub_arrays = sub_arrays;

prms.rng_seed = 42;

end % build_sim_prms_dir

%% ═══════════════════════════════════════════════════════════════════════════
