% ECE 351s - Fundamentals of Communication Systems (Fall 2026)
% (II) Digital Communication - Part A: Items 4 & 5 (Quantization and PCM)
clear; close all; clc;

fig_dir = 'digA_item4_5_figures';                 % relative folder for saved figures
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Message Signal (fs = 400 Hz = 4 x Nyquist) & True Asymmetric Range
fs = 400;                                         % Sampling frequency (Hz)
Ts = 1 / fs;
t  = 0 : Ts : 1 - Ts;                             % 1 second = 400 samples (10 periods)
m  = cos(2*pi*20*t) + 0.5*cos(2*pi*50*t);

figure('Name', 'Message Signal Plot');
plot(t, m, 'b', 'LineWidth', 1.5);
ttl = 'Message m(t)';
title(ttl);
xlabel('Time (s)');
ylabel('Amplitude (V)');
grid on;
axis([0 0.2 -1.6 1.6]);
fname = regexprep(regexprep(ttl, '[^A-Za-z0-9]+', '_'), '^_|_$', '');
print(fullfile(fig_dir, [fname '.png']), '-dpng', '-r200');

% Finely sampled copy of one period (0.1 s) to obtain true Vmax and Vmin
td   = 0 : 1/8000 : 0.1 - 1/8000;
md   = cos(2*pi*20*td) + 0.5*cos(2*pi*50*td);
Vmax = max(md);                                   % +1.5000 V
Vmin = min(md);                                   % -1.3529 V (asymmetric)
R    = Vmax - Vmin;                               %  2.8529 V
Psig = mean(m.^2);                                %  0.6250 W

%% Item 4(a) & 4(b): Mid-Tread and Mid-Rise Quantizers & Error Bound Check
n_values   = [3, 4, 6, 8];
SQNR_rise  = zeros(1, length(n_values));
SQNR_tread = zeros(1, length(n_values));

for i = 1:length(n_values)
    n = n_values(i);
    L = 2^n;

    % Mid-Tread Quantization (zero is a reconstruction level)
    Delta_tread = R / (L - 1);
    kmin  = round(Vmin / Delta_tread);
    tread = round(m / Delta_tread);
    tread(tread > kmin + L - 1) = kmin + L - 1;
    tread(tread < kmin)         = kmin;
    mq_tread  = tread * Delta_tread;
    err_tread = m - mq_tread;

    % Mid-Rise Quantization (zero is a decision threshold)
    Delta_rise = R / L;
    rise = floor((m - Vmin) / Delta_rise);
    rise(rise > L - 1) = L - 1;
    rise(rise < 0)     = 0;
    mq_rise  = Vmin + (rise + 0.5) * Delta_rise;
    err_rise = m - mq_rise;

    % Item 4(b): Verify |error| <= Delta/2
    max_err_rise  = max(abs(err_rise));
    max_err_tread = max(abs(err_tread));

    fprintf('n=%d (Rise) : Max Error = %.5f V. Bound (Delta/2 = %.5f) Valid: %d\n', ...
        n, max_err_rise, Delta_rise/2, max_err_rise <= Delta_rise/2*(1+1e-9));
    fprintf('n=%d (Tread): Max Error = %.5f V. Bound (Delta/2 = %.5f) Valid: %d\n', ...
        n, max_err_tread, Delta_tread/2, max_err_tread <= Delta_tread/2*(1+1e-9));

    % One standalone error figure per n (both types + respective +-Delta/2 bounds)
    figure('Name', ['Quantization Error (n=', num2str(n), ')']);
    h1 = plot(t, err_rise,  'r',   'LineWidth', 1.2); hold on;
    h2 = plot(t, err_tread, 'b--', 'LineWidth', 1.2);
    h3 = plot([0 0.2], [ Delta_rise/2  Delta_rise/2], 'k-',  'LineWidth', 1.5);
         plot([0 0.2], [-Delta_rise/2 -Delta_rise/2], 'k-',  'LineWidth', 1.5);
    h4 = plot([0 0.2], [ Delta_tread/2  Delta_tread/2], 'm--', 'LineWidth', 1.5);
         plot([0 0.2], [-Delta_tread/2 -Delta_tread/2], 'm--', 'LineWidth', 1.5);
    ttl = ['Quantization Error (n = ', num2str(n), ')'];
    title(ttl);
    xlabel('Time (s)'); ylabel('Error (V)');
    legend([h1 h2 h3 h4], 'Mid-Rise', 'Mid-Tread', '\pm\Delta/2 Mid-Rise', '\pm\Delta/2 Mid-Tread', ...
           'Location', 'northeast');
    dmax = max(Delta_rise, Delta_tread) / 2;
    grid on; axis([0 0.2 -1.3*dmax 1.3*dmax]);

    fname = regexprep(regexprep(ttl, '[^A-Za-z0-9]+', '_'), '^_|_$', '');
    print(fullfile(fig_dir, [fname '.png']), '-dpng', '-r200');

    SQNR_rise(i)  = 10 * log10(Psig / mean(err_rise.^2));
    SQNR_tread(i) = 10 * log10(Psig / mean(err_tread.^2));
end

%% Item 4(c): SQNR vs n (Simulation vs Theoretical Expression)
n_plot   = 3 : 0.1 : 8;
const_dB = 10 * log10(12 * Psig / R^2);                       % = -0.36 dB
theo_SQNR_rise  = 20*log10(2) .* n_plot + const_dB;           % 6.02*n - 0.36
theo_SQNR_tread = 10 * log10(Psig ./ ((R ./ (2.^n_plot - 1)).^2 / 12));

theo_vals_rise  = 20*log10(2) .* n_values + const_dB;
theo_vals_tread = 10 * log10(Psig ./ ((R ./ (2.^n_values - 1)).^2 / 12));
fprintf('\n SQNR (dB) Simulation vs Theory \n');
for i = 1:length(n_values)
    fprintf('n=%d | Rise: Sim = %.2f dB, Theo = %.2f dB | Tread: Sim = %.2f dB, Theo = %.2f dB\n', ...
        n_values(i), SQNR_rise(i), theo_vals_rise(i), SQNR_tread(i), theo_vals_tread(i));
end

figure('Name', 'Item 4c: SQNR Comparison');
plot(n_values, SQNR_rise,  'ro', 'LineWidth', 1.8, 'MarkerSize', 8); hold on;
plot(n_values, SQNR_tread, 'bs', 'LineWidth', 1.8, 'MarkerSize', 8);
plot(n_plot, theo_SQNR_rise,  'r-',  'LineWidth', 1.5);
plot(n_plot, theo_SQNR_tread, 'b--', 'LineWidth', 1.5);
ttl = 'SQNR vs Number of Bits n';
title(ttl);
xlabel('Number of Bits (n)'); ylabel('SQNR (dB)');
legend('Sim Mid-Rise', 'Sim Mid-Tread', ...
       sprintf('Theory Mid-Rise (6.02n %+.2f)', const_dB), ...
       'Theory Mid-Tread', 'Location', 'northwest');
grid on;
fname = regexprep(regexprep(ttl, '[^A-Za-z0-9]+', '_'), '^_|_$', '');
print(fullfile(fig_dir, [fname '.png']), '-dpng', '-r200');

%% Item 4(d): Minimum Bits Required for SQNR >= 40 dB
fprintf('\n Minimum Bits for 40 dB \n');
n_target = ceil((40 - const_dB) / (20 * log10(2)));           % ceil(6.70) = 7 bits
fprintf('Theoretical minimum: n = %d bits (exact = %.2f)\n', ...
        n_target, (40 - const_dB) / (20 * log10(2)));

for n_test = [n_target, n_target + 1]
    L_t = 2^n_test;

    % Mid-Rise
    Dr = R / L_t;
    r_idx = floor((m - Vmin) / Dr);
    r_idx(r_idx > L_t - 1) = L_t - 1;
    r_idx(r_idx < 0)       = 0;
    mq_r = Vmin + (r_idx + 0.5) * Dr;
    sqnr_r = 10 * log10(Psig / mean((m - mq_r).^2));

    % Mid-Tread
    Dt = R / (L_t - 1);
    km = round(Vmin / Dt);
    t_idx = round(m / Dt);
    t_idx(t_idx > km + L_t - 1) = km + L_t - 1;
    t_idx(t_idx < km)           = km;
    mq_t = t_idx * Dt;
    sqnr_t = 10 * log10(Psig / mean((m - mq_t).^2));

    fprintf('n = %d -> Sim SQNR (Rise) = %.2f dB | Sim SQNR (Tread) = %.2f dB\n', ...
            n_test, sqnr_r, sqnr_t);
end

%% Item 5: PCM Encoder, Decoder, Verification & Bandwidth Analysis
Ns    = length(m);
n_pcm = n_target;                                             % 7 bits from Item 4(d)
L_pcm = 2^n_pcm;                                              % 128 levels

% Quantize m at n_pcm = 7 bits (Mid-Tread: SQNR = 41.77 dB >= 40 dB)
Delta_pcm = R / (L_pcm - 1);
kmin_pcm  = round(Vmin / Delta_pcm);
tread_pcm = round(m / Delta_pcm);
tread_pcm(tread_pcm > kmin_pcm + L_pcm - 1) = kmin_pcm + L_pcm - 1;
tread_pcm(tread_pcm < kmin_pcm)             = kmin_pcm;

mq_pcm  = tread_pcm * Delta_pcm;
idx_pcm = tread_pcm - kmin_pcm;                               % {0, 1, ..., 127}

% Natural Binary PCM Encoder (MSB first)
codewords = zeros(Ns, n_pcm);
for b = 1:n_pcm
    weight = 2^(n_pcm - b);
    codewords(:, b) = mod(floor(idx_pcm / weight), 2);
end
pcm_bits = reshape(codewords.', 1, []);                       % 1 x 2800 serial stream

Rb = n_pcm * fs;                                              % 2800 bits/s
Tb = 1 / Rb;                                                  % 0.3571 ms
fprintf('\n PCM Encoder (n = %d bits) \n', n_pcm);
fprintf('Total encoded bits = %d bits | Bit rate Rb = %d bits/s\n', length(pcm_bits), Rb);
for k = 1:3
    fprintf('Sample %d: mq = %7.4f V | idx = %3d | bits = [%s]\n', ...
            k, mq_pcm(k), idx_pcm(k), num2str(codewords(k, :), '%d '));
end

% Plot first 3 encoded samples (21 bits)
N_show_samples = 3;
N_show_bits    = N_show_samples * n_pcm;
bits_show      = pcm_bits(1 : N_show_bits);
t_bits_ms      = (0 : N_show_bits) * Tb * 1e3;

figure('Name', 'PCM Encoded Bit Stream');
h1 = stairs(t_bits_ms, [bits_show, bits_show(end)], 'b', 'LineWidth', 1.8); hold on;
for k = 1:N_show_samples
    xb = k * n_pcm * Tb * 1e3;
    h2 = plot([xb xb], [-0.2 1.2], 'r--', 'LineWidth', 1.2);
end
ttl = 'PCM Encoded Bit Stream';
title(ttl);
xlabel('Time (ms)'); ylabel('Binary Level (0 or 1)');
legend([h1 h2], 'Encoded PCM Bits', '7-Bit Sample Boundary', 'Location', 'southeast');
grid on; axis([0 t_bits_ms(end) -0.25 1.25]);
fname = regexprep(regexprep(ttl, '[^A-Za-z0-9]+', '_'), '^_|_$', '');
print(fullfile(fig_dir, [fname '.png']), '-dpng', '-r200');

% PCM Decoder (Serial-to-Parallel + Binary-to-Decimal + Level Reconstruction)
codewords_rx = reshape(pcm_bits, n_pcm, []).';                % 400 x 7 matrix
idx_decoded  = zeros(1, Ns);
for b = 1:n_pcm
    weight = 2^(n_pcm - b);
    idx_decoded = idx_decoded + codewords_rx(:, b).' * weight;
end
mq_decoded = (idx_decoded + kmin_pcm) * Delta_pcm;

num_idx_errors = sum(idx_decoded ~= idx_pcm);
max_amp_diff   = max(abs(mq_decoded - mq_pcm));

fprintf('\n PCM Decoder Verification \n');
fprintf('Mismatched sample indices = %d (out of %d)\n', num_idx_errors, Ns);
fprintf('Max |mq_decoded - mq_pcm| = %.2e V\n', max_amp_diff);
fprintf('Identical to quantizer output? %d (1 = Yes, 0 = No)\n', max_amp_diff < 1e-12);

% Plot Quantizer Output vs Decoded Samples on the same axes
figure('Name', 'PCM Decoder Verification');
plot(td, md, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2); hold on;
stem(t(1:40), mq_pcm(1:40), 'b', 'filled', 'LineWidth', 1.2, 'MarkerSize', 5);
plot(t(1:40), mq_decoded(1:40), 'ro', 'MarkerSize', 9, 'LineWidth', 1.5);
ttl = 'Decoded Samples vs Quantizer Output';
title(ttl);
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend('Original m(t)', 'Quantizer Output m_q', 'Decoded Samples', 'Location', 'southeast');
grid on; axis([0 0.1 -1.6 1.6]);
fname = regexprep(regexprep(ttl, '[^A-Za-z0-9]+', '_'), '^_|_$', '');
print(fullfile(fig_dir, [fname '.png']), '-dpng', '-r200');

% Bit Rate and Minimum Transmission Bandwidth
W          = 50;                                              % Message bandwidth (Hz)
BW_pcm_min = Rb / 2;                                          % Minimum Nyquist PCM BW (Hz)
BW_ratio   = BW_pcm_min / W;                                  % Expansion factor (28x W)

fs_nyq     = 2 * W;                                           % Exact Nyquist rate (100 Hz)
Rb_nyq     = n_pcm * fs_nyq;
BW_nyq_min = Rb_nyq / 2;                                      % = n_pcm * W (350 Hz)

fprintf('\n Bit Rate & Bandwidth Summary \n');
fprintf('Original message bandwidth W             = %d Hz\n', W);
fprintf('At fs = %d Hz (4x Nyquist): Rb = %d bits/s | Min PCM BW = %d Hz (%.0fx W)\n', ...
        fs, Rb, BW_pcm_min, BW_ratio);
fprintf('At fs = %d Hz (1x Nyquist): Rb = %d bits/s  | Min PCM BW = %d Hz (%dx W)\n', ...
        fs_nyq, Rb_nyq, BW_nyq_min, n_pcm);
