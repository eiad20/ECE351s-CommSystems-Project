% =====================================================================
% ECE 351s - Simulation Project - Analog Part A: Items 1, 2, 3, 4 and 6
%   Item 1 : plot m(t), expected spectrum
%   Item 2 : analytical Fourier transform
%   Item 3 : numerical Fourier transform vs analytical (one graph)
%   Item 4 : 1% bandwidth of m(t)
%   Item 6 : LSB SSB (sideband filtering + phasing), comparison,
%            bandwidth, coherent demodulation
% No built-in modulation / Hilbert / filter functions are used:
% only fft, ifft, fftshift, ifftshift, sin, cos, sign, abs, real.
% Set makeFigs = 0 to run the numbers without plotting.
% =====================================================================
clear; close all; clc;
makeFigs = 1;
outDir1 = 'Analog_PA_I1to4_Figs/';  outDir6 = 'Analog_PA_I6_Figs/';
if makeFigs
  if ~exist(outDir1, 'dir'), mkdir(outDir1); end
  if ~exist(outDir6, 'dir'), mkdir(outDir6); end
end

%% ---------------- Parameters ----------------
fs = 4000;  Ts = 1/fs;          % sampling frequency  (>> 2*(fc+B))
T  = 5;     fm = 20;            % message duration and tone frequency
fc = 200;                       % carrier frequency (Item 6)
t  = 0:Ts:T;                    % 20001 points, both end points included
N  = length(t);
sincf = @(x) sin(pi*x + (x==0)*1e-30) ./ (pi*x + (x==0)*1e-30);   % own sinc

%% ---------------- Item 1 : message ----------------
m = cos(2*pi*fm*t);
m(1) = 0.5*m(1);  m(end) = 0.5*m(end);          % value 1/2 at the two jumps (t = 0 and 5)

%% ---------------- Item 3 : numerical Fourier transform ----------------
% Zero padding to 2^20 points gives df = fs/2^20 = 0.0038 Hz (the natural
% resolution 1/T = 0.2 Hz places every DFT bin on a zero of the sinc).
Nfft = 2^20;
f  = (-Nfft/2 : Nfft/2-1) * (fs/Nfft);
M  = Ts * fftshift(fft(m, Nfft));            % M(f) ~ integral of m(t)e^{-j2 pi f t}
% Item 2 : analytical transform
Ma = (T/2) * exp(-1j*pi*f*T) .* (sincf((f-fm)*T) + sincf((f+fm)*T));
err_FT = max(abs(abs(M) - abs(Ma)));  err_cplx = max(abs(M - Ma));
% naive DFT (no zero padding) for the resolution study
Mc = Ts * fftshift(fft(m(1:end-1)));  fcoarse = (-(N-1)/2 : (N-1)/2-1) * (fs/(N-1));

%% ---------------- Item 4 : 1% bandwidth ----------------
absM = abs(M);  pos = (f >= 0);
Bm = f( find(pos & (absM >= 0.01*max(absM)), 1, 'last') );   % highest f with |M|>=1%
fprintf('Item 3: complex max error = %.2e\n', err_cplx);
fprintf('Item 3: max |M|-|Ma| = %.2e   peak |M| = %.4f (theory %.4f)\n', err_FT, max(absM), T/2);
fprintf('Item 4: B_m (1%% criterion) = %.2f Hz\n', Bm);

%% ---------------- Item 6 : LSB SSB ----------------
c = cos(2*pi*fc*t);   sn = sin(2*pi*fc*t);
Nf = 2^17;                                    % padded length for filtering
ff = (-Nf/2 : Nf/2-1) * (fs/Nf);
toTime = @(X) real(ifft(ifftshift(X)));       % back to time (padded length)
cut    = @(x) x(1:N);

% DSB-SC signal (input of the sideband filter)
s_dsb = m .* c;

% (a1) sideband filtering: ideal band-pass keeping only fc-Bf <= |f| <= fc
Bf = 100;
Hlsb = double((abs(ff) >= fc-Bf) & (abs(ff) <= fc));
S_dsb = fftshift(fft(s_dsb, Nf));
s_lsb_filt = cut(toTime(S_dsb .* Hlsb));

% (a2) phasing method: own 90-degree phase shifter H(f) = -j sgn(f)
Hpi2  = -1j * sign(ff);
Mp    = fftshift(fft(m, Nf));
m_hat = cut(toTime(Mp .* Hpi2));              % Hilbert transform of m(t)
s_lsb_ph = 0.5 * (m .* c + m_hat .* sn);      % LSB: + sign ; USB would be -

% (b) comparison and bandwidth
d = s_lsb_filt - s_lsb_ph;
rms = @(x) sqrt(mean(x.^2));
in = (t > 0.1 & t < T-0.1);                      % interior: away from the two jumps
fprintf('Item 6b interior: rms diff=%.2e (%.3f%%) max diff=%.2e\n', rms(d(in)), 100*rms(d(in))/rms(s_lsb_ph(in)), max(abs(d(in))));
fprintf('Item 6b: rms(filter)=%.4f rms(phasing)=%.4f rms diff=%.2e (%.3f%%) max diff=%.2e\n', ...
        rms(s_lsb_filt), rms(s_lsb_ph), rms(d), 100*rms(d)/rms(s_lsb_ph), max(abs(d)));

% spectra on the fine grid (same method as Item 3)
spec = @(x) Ts * abs(fftshift(fft(x, Nfft)));
S_ph  = spec(s_lsb_ph);   S_fl = spec(s_lsb_filt);   S_ds = spec(s_dsb);
idxD = find(f > 0 & S_ds >= 0.01*max(S_ds));  dl = f(idxD(1));  dh = f(idxD(end));
idxP = find(f > 0 & S_ph >= 0.01*max(S_ph));  pl = f(idxP(1));  ph_ = f(idxP(end));
idxF = find(f > 0 & S_fl >= 0.01*max(S_fl));  fl = f(idxF(1));  fh = f(idxF(end));
B_dsb      = dh - dl;                % total width of the DSB-SC band
B_ssb      = fc - pl;                % sideband occupies [fc-B, fc]
fprintf('DSB-SC band: %.2f - %.2f Hz  -> B = %.2f Hz\n', dl, dh, B_dsb);
fprintf('SSB   band: %.2f - %.2f Hz (phasing) | %.2f - %.2f Hz (filter)\n', pl, ph_, fl, fh);
fprintf('SSB bandwidth (fc - f_low) = %.2f Hz ; ratio DSB/SSB = %.2f\n', B_ssb, B_dsb/B_ssb);

% (c) coherent demodulation: multiply by 4cos(wc t), low-pass at 60 Hz
Hlp = double(abs(ff) <= 60);
demod = @(s) cut(toTime(fftshift(fft(s .* (4*c), Nf)) .* Hlp));
m_rec_ph = demod(s_lsb_ph);   m_rec_fl = demod(s_lsb_filt);
m_ref    = cut(toTime(Mp .* Hlp));              % message seen through the same LPF
fprintf('Item 6c: rms err (phasing)=%.4f  (filter)=%.4f  | interior max err=%.4f | LPF-only ref err=%.2e\n', ...
        rms(m_rec_ph - m), rms(m_rec_fl - m), max(abs(m_rec_ph(in) - m(in))), rms(m_rec_ph - m_ref));
fprintf('         interior rms err (phasing)=%.4f  (filter)=%.4f\n', rms(m_rec_ph(in)-m(in)), rms(m_rec_fl(in)-m(in)));

if ~makeFigs, return; end
%% ---------------- Figures ----------------
% (the figures used in the report were drawn with the same data)
figure(1); clf;
subplot(2,1,1); plot(t, m, 'b'); grid on; xlabel('Time (s)'); ylabel('m(t)');
title('Item 1: message m(t), 0 < t < 5 s'); axis([-0.5 5.5 -1.3 1.3]);
subplot(2,1,2); plot(t, m, 'b'); grid on; xlabel('Time (s)'); ylabel('m(t)');
title('Zoom on the first 0.2 s'); axis([0 0.2 -1.3 1.3]);
print('-dpng', [outDir1 'Item1_message.png']);

figure(2); clf;
plot(f, absM, 'b', f, abs(Ma), 'r--'); grid on; xlim([-60 60]);
xlabel('Frequency (Hz)'); ylabel('|M(f)|'); legend('Numerical (FFT)', 'Analytical');
title('Items 2-3: |M(f)|, numerical vs analytical');
print('-dpng', [outDir1 'Item3_spectrum.png']);

figure(3); clf;
plot(f, absM, 'b', fcoarse, abs(Mc), 'ko'); grid on; xlim([10 30]);
xlabel('Frequency (Hz)'); ylabel('|M(f)|');
legend('Zero-padded FFT (df = 0.0038 Hz)', 'Plain DFT (df = 0.2 Hz)');
title('Item 3: effect of frequency resolution');
print('-dpng', [outDir1 'Item3_resolution.png']);

figure(4); clf;
plot(f, absM, 'b', [0 60], 0.01*max(absM)*[1 1], 'r--', Bm*[1 1], [0 3], 'k:'); grid on;
xlim([0 45]); xlabel('Frequency (Hz)'); ylabel('|M(f)|');
legend('|M(f)|', '1% of maximum', sprintf('B = %.2f Hz', Bm));
title('Item 4: 1% bandwidth of m(t)');
print('-dpng', [outDir1 'Item4_bandwidth.png']);

figure(5); clf;
subplot(2,1,1); plot(f, Ts*abs(fftshift(fft(s_dsb,Nfft))), 'b', ff, 0.4*Hlsb, 'r--'); grid on;
xlim([100 300]); xlabel('Frequency (Hz)'); ylabel('Magnitude');
legend('|S_{DSB}(f)|', 'LSB filter mask (scaled)'); title('Item 6(a): sideband filtering - input and filter');
subplot(2,1,2); plot(f, S_fl, 'b'); grid on; xlim([100 300]);
xlabel('Frequency (Hz)'); ylabel('|S_{LSB}(f)|'); title('Filter output (LSB)');
print('-dpng', [outDir6 'Item6a_SidebandFiltering.png']);

figure(6); clf;
subplot(2,1,1); plot(t, m, 'b', t, m_hat, 'r--'); grid on; xlim([0.5 0.7]);
xlabel('Time (s)'); ylabel('Amplitude'); legend('m(t)', 'm^(t) (90-degree shifted)');
title('Item 6(a): phase shifter output');
subplot(2,1,2); plot(ff, angle(Hpi2)*180/pi, 'b'); grid on; xlim([-100 100]);
xlabel('Frequency (Hz)'); ylabel('Phase of H(f) (deg)'); title('Frequency-domain 90-degree shifter');
print('-dpng', [outDir6 'Item6a_Phasing.png']);

figure(7); clf;
subplot(3,1,1); plot(t, s_lsb_filt, 'b', t, s_lsb_ph, 'r--'); grid on; xlim([2.50 2.56]);
xlabel('Time (s)'); ylabel('Amplitude'); legend('Sideband filtering', 'Phasing'); title('Item 6(b): LSB signals');
subplot(3,1,2); plot(f, S_fl, 'b', f, S_ph, 'r--'); grid on; xlim([150 210]);
xlabel('Frequency (Hz)'); ylabel('Magnitude'); legend('Sideband filtering', 'Phasing'); title('Spectra');
subplot(3,1,3); plot(t, d, 'k'); grid on; xlabel('Time (s)'); ylabel('Difference'); title('Difference of the two methods');
print('-dpng', [outDir6 'Item6b_Comparison.png']);

figure(8); clf;
plot(f, S_ds, 'b', f, S_ph, 'r'); grid on; xlim([140 260]);
xlabel('Frequency (Hz)'); ylabel('Magnitude'); legend('DSB-SC', 'LSB-SSB'); title('Item 6(b): DSB-SC vs SSB bandwidth');
print('-dpng', [outDir6 'Item6b_Bandwidth.png']);

figure(9); clf;
subplot(2,1,1); plot(t, m, 'b', t, m_rec_ph, 'r--'); grid on; xlim([0 0.3]);
xlabel('Time (s)'); ylabel('Amplitude'); legend('Original m(t)', 'Recovered'); title('Item 6(c): coherent SSB demodulation');
subplot(2,1,2); plot(t, m_rec_ph - m, 'k'); grid on; xlabel('Time (s)'); ylabel('Error'); title('Recovery error');
print('-dpng', [outDir6 'Item6c_Demodulation.png']);
