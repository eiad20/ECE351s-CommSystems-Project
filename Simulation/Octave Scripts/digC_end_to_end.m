%% ECE 351s - Digital Communication, Part C: End-to-End System Integration
    %  Chain: message -> sample -> quantize -> PCM encode -> line code ->
    %         digital modulation -> coherent detection -> decode -> reconstruct
    clear; close all; clc;
  
    % Setup figure output directory
    fig_dir = 'digC_figures';
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    set(0, 'DefaultFigurePosition', [100 100 1000 700]);
  
    %% =========================================================================
    %  1. TRANSMITTER CHAIN
    %  =========================================================================
  
    % --- Continuous-time message generation ---
    fsim = 56000;                       % Fine simulation sampling rate [Hz]
    T_total = 1.0;                      % Total duration [s]
    L_sim = fsim * T_total;
    t_sim = (0:L_sim-1) / fsim;
  
    % Analog message m(t)
    m_analog = cos(2*pi*20*t_sim) + 0.5*cos(2*pi*50*t_sim);
    W = 50;                             % Message bandwidth [Hz]                                           
                                                                                                           
    % --- Step 1: Sampling (fs = 400 Hz = 4x Nyquist) ---
    fs = 400;                           % Sampling rate [Hz]
    Ts = 1 / fs;
    N_step = round(fsim / fs);          % 140 grid points per sample                                       
    sample_indices = 1:N_step:L_sim;
    t_sampled = t_sim(sample_indices);
    m_sampled = m_analog(sample_indices);
    Ns = length(m_sampled);             % 400 samples
  
    % --- Step 2: Quantization (n = 7 bits, Mid-tread uniform quantizer) ---
    n = 7;                              % From Part A Item 4(d) (SQNR >= 40 dB)
    L_levels = 2^n;                     % 128 levels
    Vmax = 1.5;                         % Theoretical max amplitude
    Vmin = -1.3529;                     % Theoretical min amplitude
    R = Vmax - Vmin;                    % Dynamic range = 2.8529 V
    Delta = R / (L_levels - 1);         % Step size
    kmin = round(Vmin / Delta);
  
    % Mid-tread quantizer
    levels_idx = round(m_sampled / Delta);
    levels_idx(levels_idx > kmin + L_levels - 1) = kmin + L_levels - 1;
    levels_idx(levels_idx < kmin) = kmin;
    m_quantized = levels_idx * Delta;
  
    % --- Step 3: PCM Encoding ---
    % Map level indices to 0-indexed integers [0, 127]
    int_indices = levels_idx - kmin;
    tx_bits = zeros(1, Ns * n);
    for i = 1:Ns
        val = int_indices(i);
        for b = 1:n 
            tx_bits((i-1)*n + b) = bitget(val, n - b + 1);
        end
    end
    N_bits = length(tx_bits);           % 2800 bits
    Rb = N_bits / T_total;              % Bit rate = 2800 bits/s
    Tb = 1 / Rb;                        % Bit duration [s]
    spb = round(fsim * Tb);             % Samples per bit (20 samples)
  
    % --- Step 4: Line Coding (Polar NRZ) ---
    % Bit 1 -> +1 V, Bit 0 -> -1 V 
    tx_polar_nrz = kron(2*tx_bits - 1, ones(1, spb));
    tx_polar_nrz = tx_polar_nrz(1:L_sim);
  
    % --- Step 5: Digital Modulation (BPSK) --- 
    % Select carrier fc > Rb. Here fc = 14000 Hz = 5 * Rb                                                  
    fc = 14000;
    carrier_tx = cos(2*pi*fc*t_sim);
    tx_modulated = tx_polar_nrz .* carrier_tx;
  
    %% =========================================================================
    %  2. RECEIVER CHAIN
    %  =========================================================================
  
    % --- Step 6: Coherent Demodulation (BPSK Receiver) ---
    % Local oscillator (coherent: zero phase error)
    lo_rx = 2 * cos(2*pi*fc*t_sim);
    rx_mixed = tx_modulated .* lo_rx;  % Contains baseband + 2*fc component
  
    % Integrate-and-dump (matched filter) bit decision
    rx_bits = zeros(1, N_bits);
    for i = 1:N_bits
        idx_start = (i-1)*spb + 1;
        idx_end = i*spb;
        integ_val = mean(rx_mixed(idx_start:idx_end));
        rx_bits(i) = integ_val > 0;
    end
  
    % State bit errors
    num_bit_errors = sum(tx_bits ~= rx_bits);
    fprintf('Total Transmitted Bits : %d\n', N_bits); 
    fprintf('Bit Errors at Receiver : %d\n', num_bit_errors);
  
    % --- Step 7: PCM Decoding ---
    rx_int_indices = zeros(1, Ns); 
    for i = 1:Ns 
        val = 0; 
        for b = 1:n
            val = val + rx_bits((i-1)*n + b) * 2^(n - b);
        end
        rx_int_indices(i) = val;
    end
    rx_levels_idx = rx_int_indices + kmin;
    m_decoded = rx_levels_idx * Delta;
  
    % --- Step 8: Reconstruction Filter (Ideal Low-Pass Filter) ---
    % Create sampled pulse train on fine simulation grid
    m_impulse = zeros(1, L_sim); 
    m_impulse(sample_indices) = m_decoded;
  
    % Ideal LPF with cutoff fc_lpf = 50 Hz, gain = N_step 
    f_axis = (-L_sim/2 : L_sim/2 - 1) * (fsim / L_sim); 
    H_lpf = double(abs(f_axis) <= W) * N_step;
    M_impulse_fft = fftshift(fft(m_impulse));
    m_reconstructed = real(ifft(ifftshift(M_impulse_fft .* H_lpf)));
  
    %% =========================================================================
    %  3. ERROR MEASUREMENT & ANALYSIS 
    %  =========================================================================
    % Reconstruction error vs original analog message
    % Ignore 5% edge transients caused by finite-duration ideal LPF sinc tails 
    edge_margin = round(0.05 * L_sim); 
    eval_indices = edge_margin : (L_sim - edge_margin);
  
    recon_error = m_analog(eval_indices) - m_reconstructed(eval_indices); 
    recon_mse = mean(recon_error.^2); 
    quant_mse = mean((m_sampled - m_quantized).^2); 
  
    fprintf('Quantization MSE     : %.6e V^2\n', quant_mse);
    fprintf('Reconstruction MSE   : %.6e V^2\n', recon_mse);
    fprintf('Reconstruction SQNR  : %.2f dB\n', 10*log10(mean(m_analog(eval_indices).^2)/recon_mse));
  
    %% =========================================================================
    %  4. PLOTTING AND EXPORTING FIGURES
    %  =========================================================================
  
    % --- Figure 1: Transmitter Block Signals (Time domain zoom) ---
    figure(1); 
    tw_zoom = t_sim <= 0.02;  % First 20 ms
    tw_samp = t_sampled <= 0.02; 
    subplot(4,1,1);
    plot(t_sim(tw_zoom)*1000, m_analog(tw_zoom), 'b', 'LineWidth', 1.2); hold on;
    stem(t_sampled(tw_samp)*1000, m_sampled(tw_samp), 'r', 'filled');
    grid on; ylabel('Amp (V)'); title('1. Original Analog Message and Sampled Points');
    legend('m(t)', 'Samples (fs=400Hz)', 'Location', 'eastoutside');
  
    subplot(4,1,2); 
    stem(t_sampled(tw_samp)*1000, m_quantized(tw_samp), 'm', 'filled');
    grid on; ylabel('Amp (V)'); title('2. Quantized Samples (Mid-tread, n=7 bits)');
  
    subplot(4,1,3); 
    plot(t_sim(tw_zoom)*1000, tx_polar_nrz(tw_zoom), 'k', 'LineWidth', 1.2); 
    grid on; ylabel('Amp (V)'); ylim([-1.5 1.5]); title('3. Polar NRZ Line Code');
  
    subplot(4,1,4); 
    plot(t_sim(tw_zoom)*1000, tx_modulated(tw_zoom), 'Color', [0 0.5 0], 'LineWidth', 0.8);
    grid on; xlabel('Time (ms)'); ylabel('Amp (V)'); title('4. BPSK Modulated Waveform (fc = 14 kHz)');
    print('-dpng', '-r200', fullfile(fig_dir, '01_tx_blocks_waveform.png'));
  
    % --- Figure 2: Receiver Detection & Bit Verification ---
    figure(2);
    subplot(2,1,1); 
    plot(t_sim(tw_zoom)*1000, rx_mixed(tw_zoom), 'Color', [0.6 0.6 0.6]); hold on;
    plot(t_sim(tw_zoom)*1000, tx_polar_nrz(tw_zoom), 'b', 'LineWidth', 1.5);
    grid on; xlabel('Time (ms)'); ylabel('Amp (V)');
    title('Coherent Detection: Mixer Output vs Transmitted Polar NRZ');
    legend('Mixer Output (Baseband + 2fc)', 'Transmitted Baseband NRZ', 'Location', 'northeast');
  
    subplot(2,1,2);
    n_show_bits = 40; 
    stairs(0:n_show_bits-1, tx_bits(1:n_show_bits), 'b', 'LineWidth', 2); hold on;
    stairs(0:n_show_bits-1, rx_bits(1:n_show_bits), 'r--', 'LineWidth', 1.5);
    grid on; ylim([-0.25 1.25]); xlabel('Bit Index'); ylabel('Binary Value'); 
    title(sprintf('Transmitted vs Recovered Bits (Total Bit Errors: %d / %d)', num_bit_errors, N_bits));
    legend('Transmitted Bits', 'Received Bits', 'Location', 'southeast');
    print('-dpng', '-r200', fullfile(fig_dir, '02_rx_bit_recovery.png')); 
  
    % --- Figure 3: Reconstructed Analog Message vs Original ---                                           
    figure(3);
    subplot(2,1,1);
    tw_plot = t_sim >= 0.05 & t_sim <= 0.25;
    plot(t_sim(tw_plot), m_analog(tw_plot), 'b', 'LineWidth', 1.8); hold on;
    plot(t_sim(tw_plot), m_reconstructed(tw_plot), 'r--', 'LineWidth', 1.5);
    grid on; xlabel('Time (s)'); ylabel('Amplitude (V)');
    title('End-to-End Reconstructed Signal vs Original Message');
    legend('Original m(t)', 'Reconstructed m_r(t)', 'Location', 'southeast');
  
    subplot(2,1,2);
    plot(t_sim(eval_indices), recon_error, 'm', 'LineWidth', 1.0);
    grid on; xlabel('Time (s)'); ylabel('Error (V)');
    title(sprintf('Final Reconstruction Error e(t) = m(t) - m_r(t) (MSE = %.2e V^2)', recon_mse));         
    print('-dpng', '-r200', fullfile(fig_dir, '03_reconstruction_and_error.png'));
  
    % --- Figure 4: Spectra Across the Chain ---
    figure(4);
    f_plot = (-fsim/2:fsim/2-1);
    subplot(3,1,1);
    plot(f_axis/1000, abs(fftshift(fft(m_analog))) / L_sim, 'b');
    grid on; xlim([-1 1]); xlabel('Frequency (kHz)'); ylabel('Magnitude');
    title('Spectrum of Original Message m(t)');
  
    subplot(3,1,2);
    plot(f_axis/1000, abs(fftshift(fft(tx_modulated))) / L_sim, 'Color', [0 0.5 0]);
    grid on; xlim([-20 20]); xlabel('Frequency (kHz)'); ylabel('Magnitude');
    title('Transmitted BPSK Modulated Spectrum centered at fc = 14 kHz');
  
    subplot(3,1,3); 
    plot(f_axis/1000, abs(fftshift(fft(m_reconstructed))) / L_sim, 'r');
    grid on; xlim([-1 1]); xlabel('Frequency (kHz)'); ylabel('Magnitude');
    title('Spectrum of Reconstructed Signal after Ideal LPF');
    print('-dpng', '-r200', fullfile(fig_dir, '04_system_spectra.png'));
  
    fprintf('\nAll figures successfully saved in %s/\n', fig_dir);
