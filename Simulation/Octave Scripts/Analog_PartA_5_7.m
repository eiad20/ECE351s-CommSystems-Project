% =========================================================================
% ECE 351s - Fundamentals of Communication Systems (Fall 2026)
% Simulation Project: Analog Part A - Items 5 & 7
% Author: Student Contribution (Analog Part A Items 5 & 7)
% =========================================================================
clear all; close all; clc;

%% -------------------------------------------------------------------------
% GLOBAL PARAMETERS & MESSAGE SIGNAL GENERATION
% -------------------------------------------------------------------------
fm = 20;               % Message frequency (Hz)
fc = 200;              % Carrier frequency (Hz)
T_total = 5;           % Signal duration (seconds)
fs = 4000;             % Sampling frequency (Hz) - well above Nyquist
dt = 1/fs;             % Sampling interval (seconds)
t = 0:dt:(T_total - dt);% Time vector (0 to 5 seconds)
N = length(t);         % Total number of samples

% Defined message signal m(t) = cos(2*pi*20*t) for 0 < t < 5 s
m_t = cos(2 * pi * fm * t);

% Frequency vector for FFT plotting (in Hz, centered at 0 Hz)
f = linspace(-fs/2, fs/2, N);

%% =========================================================================
% ITEM 5: DOUBLE SIDEBAND - SUPPRESSED CARRIER (DSB-SC) MODULATION
% =========================================================================
fprintf('--- Running Item 5: DSB-SC Modulation & Coherent Detection ---\n');

% Item 5.a: DSB-SC Modulation
c_t = cos(2 * pi * fc * t);          % Carrier signal
s_dsbsc = m_t .* c_t;                % DSB-SC Modulated Signal

% Compute Spectrum using basic FFT math
S_dsbsc_fft = fftshift(fft(s_dsbsc)) / N;
S_dsbsc_mag = abs(S_dsbsc_fft);

% Estimate Bandwidth (1% max threshold criterion)
max_mag_dsb = max(S_dsbsc_mag);
idx_dsb = find(S_dsbsc_mag >= 0.01 * max_mag_dsb);
f_dsb_active = abs(f(idx_dsb));
bw_dsbsc_sim = 2 * (max(f_dsb_active) - fc); % BW = 2 * fm
fprintf('DSB-SC Estimated Bandwidth: %.2f Hz (Theoretical: %.2f Hz)\n', ...
    bw_dsbsc_sim, 2*fm);

% Plotting Item 5.a Waveform & Spectrum
figure('Name', 'Item 5a: DSB-SC Signal & Spectrum');
subplot(2,1,1);
plot(t, s_dsbsc, 'b', 'LineWidth', 1); hold on;
plot(t, m_t, 'r--', 'LineWidth', 1.2);
plot(t, -m_t, 'r--', 'LineWidth', 1.2);
title('Item 5(a): DSB-SC Modulated Waveform s_{DSB}(t) with Message Envelope');
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend('DSB-SC Signal s(t)', 'Envelope \pm m(t)'); grid on; xlim([0 0.25]);

subplot(2,1,2);
plot(f, S_dsbsc_mag, 'b', 'LineWidth', 1.2);
title('Item 5(a): Magnitude Spectrum |S_{DSB}(f)|');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
grid on; xlim([-300 300]);

% -------------------------------------------------------------------------
% Item 5.b: Coherent Demodulation (Ideal Phase: phi = 0)
% -------------------------------------------------------------------------
% Multiply received signal by local oscillator 2*cos(w_c*t)
v_mult_ideal = s_dsbsc .* (2 * cos(2 * pi * fc * t));

% Frequency-domain Ideal Low Pass Filter (Cutoff at W = 40 Hz)
V_ideal_fft = fftshift(fft(v_mult_ideal)) / N;
LPF_mask = abs(f) <= (fm + 10); % Pass frequencies up to 30 Hz
m_rec_ideal = real(ifft(ifftshift(V_ideal_fft .* LPF_mask))) * N;

figure('Name', 'Item 5b: Coherent Detection (Ideal Phase)');
plot(t, m_t, 'b', 'LineWidth', 1.5); hold on;
plot(t, m_rec_ideal, 'r--', 'LineWidth', 1.2);
title('Item 5(b): Coherent Detector Output (\phi = 0^\circ) vs Original Message');
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend('Original Message m(t)', 'Recovered Message m_{rec}(t)');
grid on; xlim([0 0.25]);

% -------------------------------------------------------------------------
% Item 5.c: Phase Error Analysis (phi = 30 deg, 60 deg, 90 deg)
% -------------------------------------------------------------------------
phases_deg = [30, 60, 90];
m_rec_phase = zeros(length(phases_deg), N);

figure('Name', 'Item 5c: Coherent Detection with Phase Errors');
plot(t, m_t, 'k', 'LineWidth', 1.5); hold on;

for k = 1:length(phases_deg)
    phi_rad = deg2rad(phases_deg(k));
    v_mult_phase = s_dsbsc .* (2 * cos(2 * pi * fc * t + phi_rad));
    V_phase_fft = fftshift(fft(v_mult_phase)) / N;
    m_rec_phase(k, :) = real(ifft(ifftshift(V_phase_fft .* LPF_mask))) * N;

    plot(t, m_rec_phase(k, :), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Demodulated (\\phi = %d^\\circ)', phases_deg(k)));
end

title('Item 5(c): Effect of Local Oscillator Phase Error on Coherent Detection');
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend('Original m(t)', '\phi = 30^\circ', '\phi = 60^\circ', '\phi = 90^\circ (Quadrature Null)');
grid on; xlim([0 0.25]);

%% =========================================================================
% ITEM 7: CONVENTIONAL AMPLITUDE MODULATION (AM) & ENVELOPE DETECTOR
% =========================================================================
fprintf('\n--- Running Item 7: Conventional AM & Envelope Detector ---\n');

% -------------------------------------------------------------------------
% Item 7.a & 7.b: Normal Modulation (ka = 0.4, mu = 0.4)
% -------------------------------------------------------------------------
ka1 = 0.4;
s_am1 = (1 + ka1 * m_t) .* cos(2 * pi * fc * t);

% Compute Spectrum
S_am1_fft = fftshift(fft(s_am1)) / N;
S_am1_mag = abs(S_am1_fft);

% Envelope Detector Algorithm (Diode + RC Circuit Simulation)
% Discrete update: V_out[n] = max(s_am[n], V_out[n-1] * exp(-dt/RC))
RC_opt = 0.015; % Chosen such that 1/fc (0.005s) << RC << 1/fm (0.05s) and RC <= sqrt(1-mu^2)/(2*pi*fm*mu)
v_env1 = zeros(1, N);
v_out = 0;
decay_factor_opt = exp(-dt / RC_opt);

for n = 1:N
    v_in = s_am1(n);
    if v_in > v_out
        v_out = v_in; % Diode conducts (charge)
    else
        v_out = v_out * decay_factor_opt; % Diode cuts off (discharge)
    end
    v_env1(n) = v_out;
end

% Remove DC Component to recover AC Message
v_env1_ac = v_env1 - mean(v_env1);
v_env1_ac_scaled = v_env1_ac / ka1; % Scale back amplitude

figure('Name', 'Item 7a & 7b: Conventional AM (ka = 0.4)');
subplot(2,1,1);
plot(f, S_am1_mag, 'r', 'LineWidth', 1.2);
title('Item 7(a): Spectrum of Conventional AM Signal (k_a = 0.4)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
grid on; xlim([-300 300]);

subplot(2,1,2);
plot(t, s_am1, 'g', 'LineWidth', 0.8); hold on;
plot(t, v_env1, 'k', 'LineWidth', 1.5);
plot(t, (1 + ka1 * m_t), 'r--', 'LineWidth', 1.2);
title(sprintf('Item 7(b): AM Signal (k_a = 0.4) & Envelope Detector Output (RC = %.3f s)', RC_opt));
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend('AM Wave s_{AM}(t)', 'Envelope Detector Output', 'True Envelope (1 + k_a m(t))');
grid on; xlim([0 0.25]);

% -------------------------------------------------------------------------
% RC Relaxation Time Comparison (Too Small, Suitable, Too Large)
% -------------------------------------------------------------------------
RC_small = 0.001;  % RC << 1/fc (Excessive ripple)
RC_large = 0.100;  % RC > sqrt(1-mu^2)/(2*pi*fm*mu) (Diagonal clipping)

v_env_small = zeros(1, N); v_env_large = zeros(1, N);
v_out_s = 0; v_out_l = 0;
decay_small = exp(-dt / RC_small); decay_large = exp(-dt / RC_large);

for n = 1:N
    % Small RC
    if s_am1(n) > v_out_s, v_out_s = s_am1(n); else, v_out_s = v_out_s * decay_small; end
    v_env_small(n) = v_out_s;
    % Large RC
    if s_am1(n) > v_out_l, v_out_l = s_am1(n); else, v_out_l = v_out_l * decay_large; end
    v_env_large(n) = v_out_l;
end

figure('Name', 'Item 7b: Effect of RC Time Constant on Demodulation');
plot(t, s_am1, 'Color', [0.7 0.7 0.7]); hold on;
plot(t, v_env_small, 'b', 'LineWidth', 1.2);
plot(t, v_env1, 'g', 'LineWidth', 1.5);
plot(t, v_env_large, 'r', 'LineWidth', 1.5);
title('Item 7(b): Envelope Detector Output for Different RC Values (k_a = 0.4)');
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend('AM Signal', 'Too Small RC (0.001 s - High Ripple)', ...
    'Optimal RC (0.015 s - Good Tracking)', 'Too Large RC (0.100 s - Diagonal Clipping)');
grid on; xlim([0 0.25]);

% -------------------------------------------------------------------------
% Item 7.c: Over-Modulation (ka = 1.5, mu = 1.5)
% -------------------------------------------------------------------------
ka2 = 1.5;
s_am2 = (1 + ka2 * m_t) .* cos(2 * pi * fc * t);

% Compute Spectrum for ka = 1.5
S_am2_fft = fftshift(fft(s_am2)) / N;
S_am2_mag = abs(S_am2_fft);

% Demodulate Over-modulated Signal with Envelope Detector (Optimal RC)
v_env2 = zeros(1, N);
v_out = 0;
for n = 1:N
    v_in = s_am2(n);
    if v_in > v_out
        v_out = v_in;
    else
        v_out = v_out * decay_factor_opt;
    end
    v_env2(n) = v_out;
end

figure('Name', 'Item 7c: Over-modulation Analysis (ka = 1.5)');
subplot(2,1,1);
plot(f, S_am2_mag, 'r', 'LineWidth', 1.2);
title('Item 7(c): Spectrum of Over-modulated AM Signal (k_a = 1.5)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
grid on; xlim([-300 300]);

subplot(2,1,2);
plot(t, s_am2, 'g', 'LineWidth', 0.8); hold on;
plot(t, v_env2, 'k', 'LineWidth', 1.5);
plot(t, abs(1 + ka2 * m_t), 'r--', 'LineWidth', 1.2);
title('Item 7(c): Over-modulation Waveform (k_a = 1.5) & Envelope Detector Distortion');
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend('Over-modulated AM Wave', 'Envelope Detector Output', 'Rectified True Envelope |1 + k_a m(t)|');
grid on; xlim([0 0.25]);

fprintf('--- Script Execution Completed Successfully ---\n');
