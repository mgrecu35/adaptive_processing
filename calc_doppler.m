function X_dop=calc_doppler(X,prms)
    N  = prms.N_pulses;
    win = hann(N);
    win = win / mean(win);           % amplitude normalisation

% Windowed DFT along pulse dimension [N_ch x N_fft]
    X_dop = zeros(prms.N_ch, N);
    for m = 1:prms.N_ch
        X_dop(m,:) = fft(X(m,:) .* win.', N);
    end
end