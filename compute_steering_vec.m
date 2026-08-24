function s = compute_steering_vec(rx_pos, az_rad, el_rad, lambda)
% rx_pos : [3 x M] subarray phase centres (from meta.rx_pos)
% az_rad : azimuth of look direction [rad]
% el_rad : elevation of look direction [rad]
% lambda : wavelength [m]

u = [cos(el_rad)*cos(az_rad);   % boresight component
     cos(el_rad)*sin(az_rad);   % lateral
     sin(el_rad)];              % vertical

phi = (2*pi/lambda) * (u.' * rx_pos);  % [1 x M]
s   = exp(1j * phi).';                 % [M x 1]
end
