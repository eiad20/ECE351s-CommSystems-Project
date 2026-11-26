%% ECE 351s - Digital Communication, Part A, Items 1-3
%  Sampling (Nyquist / 4x Nyquist), aliasing, natural & flat-top sampling
%  Works in Octave and MATLAB. No toolbox communication functions are used:
%  every block (sampler, ideal LPF, pulse train, equalizer) is built here
%  from basic theory.  Run the whole file; figures are also saved as PNGs
%  in a "figures" folder next to this script (relative path, no editing).

clear; close all; clc;
saveFigs = true;
set(0,'DefaultFigurePosition',[100 100 1150 780]);      % bigger figures so nothing is cramped
set(0,'DefaultFigurePaperPositionMode','auto');
% legends go OUTSIDE the axes (right side) so they never cover the curves
leg = @(varargin) legend(varargin{:},'Location','eastoutside');

%% ---------------------------------------------------------------------
%  0. SETUP: fake "continuous time" with a very fine grid
%  ---------------------------------------------------------------------
%  fsim = 16000 Hz and T = 1 s are chosen so that
%   * 16000 is divisible by 80, 100 and 400  -> every sampling rate is an
%     exact integer number of grid points (N = fsim/fs)
%   * 25 % of the 400 Hz period (2.5 ms) is exactly 10 grid points
%   * 1 s holds a whole number of cycles of 20 Hz and 50 Hz, so the FFT has
%     no leakage and 1 Hz resolution: every spectral line is a single bin
fsim = 16000;                         % fine grid rate [Hz]
T    = 1;                             % observation window [s]
L    = fsim*T;                        % number of grid points
t    = (0:L-1)/fsim;                  % time axis [s]
f    = -fsim/2 : 1/T : fsim/2-1/T;    % frequency axis [Hz], matches fftshift
W    = 50;                            % message bandwidth [Hz]
fNyq = 2*W;                           % Nyquist rate = 100 Hz

m = cos(2*pi*20*t) + 0.5*cos(2*pi*50*t);     % the message

% small helpers (anonymous functions so the file is MATLAB/Octave portable)
mysinc = @(x) (x==0) + (x~=0).*sin(pi*x)./(pi*x + (x==0));   % sin(pi x)/(pi x)
spec   = @(x) fftshift(fft(x))/length(x);    % centred spectrum, line height = Fourier coefficient
applyH = @(x,H) real(ifft(ifftshift(fftshift(fft(x)).*H)));  % filter in frequency domain
savefig_ = @(name) print('-dpng','-r110',fullfile('figures',[name '.png']));
if saveFigs && ~exist('figures','dir'), mkdir('figures'); end

% reusable measurement: height of the spectral line at frequency f0 [Hz]
lineAt = @(S,f0) abs(S(f==f0));

% ---- message and its spectrum ----
figure('Name','Message');
subplot(2,1,1); plot(t,m,'b','LineWidth',1.2); grid on; xlim([0 0.2]);
xlabel('Time (s)'); ylabel('m(t) (V)'); title('Message m(t) = cos(2\pi20t) + 0.5cos(2\pi50t)');
subplot(2,1,2); plot(f,abs(spec(m)),'b','LineWidth',1.2); grid on; xlim([-100 100]); ylim([0 0.6]);
xlabel('Frequency (Hz)'); ylabel('|M(f)|'); title('Message spectrum: lines at \pm20 Hz (0.5) and \pm50 Hz (0.25)');
if saveFigs, savefig_('00_message'); end

%% ---------------------------------------------------------------------
%  ITEM 1: sampling at Nyquist (100 Hz) and 4x Nyquist (400 Hz)
%  ---------------------------------------------------------------------
%  Ideal sampling = multiplying by an impulse train. On the fine grid we
%  keep every N-th point and set the others to zero ("zero-stuffing").
%  Theory:  Ms(f) = fs * sum_k M(f - k fs)   -> copies of M(f) at k*fs.
%  The zero-stuffed signal has spectrum (1/N)*sum_k M(f-k fs), so we
%  multiply its spectrum by N when plotting (copies then have the same
%  height as M(f)) and use gain N in the reconstruction LPF.
rates1 = [fNyq, 4*fNyq];
res1 = struct();
for i = 1:2
    fs = rates1(i);  N = fsim/fs;
    ms = zeros(size(m));  ms(1:N:end) = m(1:N:end);       % sampled signal
    Ms = N*spec(ms);                                       % spectrum with copies

    % ---- own ideal low-pass filter (frequency-domain brick wall, gain N)
    fc = fs/2;                                             % cutoff = fs/2
    if fs == fNyq
        H = N*((abs(f) <  fc) + 0.5*(abs(f) == fc));       % value 1/2 AT the edge (ideal-filter convention)
    else
        H = N*(abs(f) <= fc);
    end
    mr = applyH(ms,H);

    % ---- second, independent LPF: windowed-sinc FIR with conv (check)
    M_ = 400; n_ = -M_:M_;                                 % 801 taps
    h  = (2*fc/fsim)*mysinc(2*fc*n_/fsim) .* (0.54+0.46*cos(pi*n_/M_));  % Hamming window by hand
    h  = h/sum(h);                                         % unity DC gain
    mr_fir = N*conv(ms,h,'same');
    mid = 2*M_+1 : L-2*M_;                                 % ignore FIR edge transients
    res1(i).fs = fs;  res1(i).maxErr = max(abs(mr-m));
    res1(i).maxErrFIR = max(abs(mr_fir(mid)-m(mid)));

    figure('Name',sprintf('Item1 fs=%d',fs));
    tw = t<=0.1;
    subplot(3,1,1); plot(t(tw),m(tw),'b'); hold on;
    idx = find(ms~=0 & tw); stem(t(idx),ms(idx),'r','filled','MarkerSize',3);
    grid on; xlabel('Time (s)'); ylabel('Amplitude (V)'); leg('m(t)','samples');
    if fs==fNyq, lbl='Nyquist rate'; else, lbl='4 x Nyquist'; end
    title(sprintf('Sampled signal, f_s = %d Hz (%s)',fs,lbl));
    subplot(3,1,2); plot(f,abs(Ms),'r'); grid on;
    xlim([-3.5*fs 3.5*fs]); ylim([0 0.65]); xlabel('Frequency (Hz)'); ylabel('|M_s(f)|');
    title(sprintf('Spectrum of sampled signal: copies of M(f) every %d Hz',fs));
    subplot(3,1,3); plot(t(tw),m(tw),'b','LineWidth',1.4); hold on;
    plot(t(tw),mr(tw),'r--','LineWidth',1.4); plot(t(tw),mr_fir(tw),'g:','LineWidth',1.4);
    grid on; xlabel('Time (s)'); ylabel('Amplitude (V)');
    leg('original m(t)','ideal LPF (freq. domain)','windowed-sinc FIR (approximate)');
    title({sprintf('Reconstruction, f_s = %d Hz',fs), sprintf('max error: ideal LPF = %.1e, FIR (approximate) = %.1e',res1(i).maxErr,res1(i).maxErrFIR)});
    if saveFigs, savefig_(sprintf('01_item1_fs%d',fs)); end
end

% ---- The 50 Hz component EXACTLY at the Nyquist rate -------------------
%  Samples of 0.5cos(2*pi*50*t) at fs=100 are 0.5(-1)^n.  The 50 Hz line of
%  the base-band copy and of the neighbouring copy land on the SAME bin, and
%  the LPF edge sits exactly on that bin.  The result depends on the value
%  the "ideal" filter takes AT its edge: 1, 1/2, or 0.
fs = fNyq; N = fsim/fs;
ms = zeros(size(m)); ms(1:N:end) = m(1:N:end);
edgeGain = [1 0.5 0];
figure('Name','Item1 edge');
tw = t<=0.1; col = {'r','g','m'};
subplot(2,1,1); plot(t(tw),m(tw),'b','LineWidth',1.6); hold on;
for k = 1:3
    H  = N*((abs(f) < fs/2) + edgeGain(k)*(abs(f) == fs/2));
    mk = applyH(ms,H);
    A50(k) = 2*lineAt(spec(mk),50);                        % recovered amplitude of the 50 Hz cosine
    plot(t(tw),mk(tw),[col{k} '--'],'LineWidth',1.2);
end
grid on; xlabel('Time (s)'); ylabel('Amplitude (V)');
leg('original','edge gain = 1','edge gain = 1/2','edge gain = 0');
title('f_s = 100 Hz: result depends on LPF gain at exactly 50 Hz');

% sine of the same frequency: all its samples are (numerically) zero
msin   = sin(2*pi*50*t);
subplot(2,1,2); plot(t(tw),msin(tw),'b','LineWidth',1.4); hold on;
idxs = find(mod(0:L-1,N)==0 & t<=0.1);
stem(t(idxs),msin(idxs),'r','filled','MarkerSize',3);
grid on; xlabel('Time (s)'); ylabel('Amplitude (V)');
leg('sin(2\pi50t)','samples at 100 Hz');
title('sin(2\pi50t) sampled at 100 Hz: every sample is 0');
if saveFigs, savefig_('02_item1_nyquist_edge'); end
sineMaxSample = max(abs(msin(1:N:end)));

%% ---------------------------------------------------------------------
%  ITEM 2: sampling BELOW the Nyquist rate (fs = 80 Hz < 100 Hz)
%  ---------------------------------------------------------------------
%  80 Hz is still > 2*20, so the 20 Hz tone is sampled legally, but the
%  50 Hz tone is not:  a tone f0 reappears at |f0 - k fs|  ->  |50-80| = 30 Hz.
fs = 80; N = fsim/fs;
ms = zeros(size(m)); ms(1:N:end) = m(1:N:end);
Ms = N*spec(ms);
fc = 40;                                       % best LPF we can use: fs/2
mr = applyH(ms, N*(abs(f)<=fc));
Mr = spec(mr);
peaks = f(abs(Mr)>0.01 & f>0);                 % frequencies present after reconstruction
mAlias = cos(2*pi*20*t) + 0.5*cos(2*pi*30*t);  % what theory predicts we get
errAlias = max(abs(mr-mAlias));                % should be ~0: output IS the aliased signal
errAliasVsOrig = sqrt(mean((mr-m).^2));

figure('Name','Item2 aliasing');
subplot(3,1,1); plot(f,abs(Ms),'r'); hold on; grid on; xlim([-260 260]); ylim([0 0.65]);
yl = ylim; plot([fc fc],yl,'k:'); plot([-fc -fc],yl,'k:');
text(31,0.38,'\leftarrow alias','Color','m');
xlabel('Frequency (Hz)'); ylabel('|M_s(f)|');
title('f_s = 80 Hz < Nyquist: copies overlap, 50 Hz folds to 30 Hz');
subplot(3,1,2); plot(f,abs(Mr),'b'); grid on; xlim([-100 100]);
xlabel('Frequency (Hz)'); ylabel('|M_r(f)|');
title('Spectrum after ideal LPF (fc = 40 Hz): lines at 20 Hz and 30 Hz, none at 50 Hz');
subplot(3,1,3); tw = t<=0.2;
plot(t(tw),m(tw),'b','LineWidth',1.4); hold on; plot(t(tw),mr(tw),'r--','LineWidth',1.4);
grid on; xlabel('Time (s)'); ylabel('Amplitude (V)'); leg('original m(t)','reconstructed');
title(sprintf('Cannot be recovered: RMS error = %.3f V',errAliasVsOrig));
if saveFigs, savefig_('03_item2_aliasing'); end

%% ---------------------------------------------------------------------
%  ITEM 3: natural vs flat-top sampling, 25 % duty cycle, fs = 400 Hz
%  ---------------------------------------------------------------------
fs = 4*fNyq; N = fsim/fs; Ts = 1/fs;           % Ts = 2.5 ms, N = 40 grid points
d  = 0.25;  w = round(d*N);  tau = w/fsim;     % pulse = 10 points = 0.625 ms

% -- natural sampling: m(t) * p(t), p(t) = periodic pulse train
p      = double(mod(0:L-1,N) < w);
m_nat  = m .* p;
% -- flat-top: read m ONLY at the sampling instant, hold it for tau
m_flat = zeros(size(m));
for n = 0:fs*T-1
    m_flat(n*N+1 : n*N+w) = m(n*N+1);
end
S_nat  = spec(m_nat);   S_flat = spec(m_flat);

% -- recovery
H_lp    = (abs(f) <= fs/2);                         % ideal LPF, gain 1
mr_nat  = applyH(m_nat, H_lp/d);                    % natural: LPF + gain 1/d
mr_flat = applyH(m_flat,H_lp/d);                    % flat-top: same filter -> distorted
% equalizer: undo sinc(f*tau) magnitude (droop) and the tau/2 hold delay
delay   = (w-1)/(2*fsim);                           % centre of a w-point hold on the grid
Heq     = zeros(size(f)); in = abs(f)<=W;
Heq(in) = exp(1j*2*pi*f(in)*delay) ./ (d*mysinc(f(in)*tau));
mr_eq   = applyH(m_flat,Heq);
% flat-top without equalizer: remove only the delay so the comparison is fair
Hd      = zeros(size(f)); Hd(in) = exp(1j*2*pi*f(in)*delay)/d;
mr_flat_d = applyH(m_flat,Hd);

% -- aperture-effect numbers at the band edge
droop50   = mysinc(W*tau);        droop50dB = 20*log10(droop50);
droop20   = mysinc(20*tau);
boost50dB = -droop50dB;
droopFull = mysinc(W*Ts);         droopFulldB = 20*log10(droopFull);  % 100 % duty, for context
ratio20 = lineAt(S_flat,20)/lineAt(S_nat,20);
ratio50 = lineAt(S_flat,50)/lineAt(S_nat,50);

figure('Name','Item3 time');
tw = t<=0.008;
plot(t(tw),m(tw),'b','LineWidth',1.2); hold on;
stairs(t(tw),m_nat(tw),'r','LineWidth',1.8);
stairs(t(tw),m_flat(tw),'g--','LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Amplitude (V)');
leg('m(t)','natural sampled','flat-top sampled');
title({'Natural vs flat-top sampling (25 % duty, f_s = 400 Hz)','natural tops follow m(t), flat-top tops are constant'});
if saveFigs, savefig_('04_item3_time'); end

figure('Name','Item3 spectra');
env = 0.5*d*abs(mysinc(f*tau));                     % sinc envelope for the 20 Hz lines
plot(f,abs(S_nat),'r','LineWidth',1.8); hold on;
plot(f,abs(S_flat),'g--','LineWidth',1.3);
plot(f,env,'k:','LineWidth',1.2); grid on; xlim([-2000 2000]); ylim([0 0.15]);
xlabel('Frequency (Hz)'); ylabel('|S(f)|');
leg('natural','flat-top','envelope d\cdot sinc(f\tau)\cdot0.5');
title({'Natural vs flat-top spectra (f_s = 400 Hz)','natural copies are constant-scaled, flat-top copies are also multiplied by sinc(f\tau)'});
if saveFigs, savefig_('05_item3_spectra'); end

figure('Name','Item3 recovery','Position',[50 30 1250 900]);   % wider/taller so titles, legends and labels fit
tw = t<=0.1;
subplot('Position',[0.10 0.72 0.60 0.20]); plot(t(tw),m(tw),'b','LineWidth',1.6); hold on;
plot(t(tw),mr_nat(tw),'r--','LineWidth',1.2); plot(t(tw),mr_flat_d(tw),'g--','LineWidth',1.2);
plot(t(tw),mr_eq(tw),'m:','LineWidth',1.6);
grid on; xlabel('Time (s)'); ylabel('Amplitude (V)');
leg('original','natural + LPF','flat-top + LPF (no equalizer)','flat-top + LPF + equalizer');
title({'Recovery from natural and flat-top samples','(curves overlap: the effect is small)'});
% error panel: difference from the original, so the equalizer's benefit is visible
subplot('Position',[0.10 0.40 0.60 0.20]);
plot(t(tw),mr_nat(tw)-m(tw),'r','LineWidth',1.2); hold on;
plot(t(tw),mr_flat_d(tw)-m(tw),'g','LineWidth',1.2);
plot(t(tw),mr_eq(tw)-m(tw),'m','LineWidth',1.2);
grid on; xlabel('Time (s)'); ylabel('Error (V)');
leg('natural + LPF','flat-top, no equalizer','flat-top + equalizer');
title('Recovery error (reconstructed - original)');
subplot('Position',[0.10 0.08 0.60 0.20]); ff = 0:1:W;
plot(ff,20*log10(mysinc(ff*tau)),'r','LineWidth',1.4); hold on;
plot(ff,-20*log10(mysinc(ff*tau)),'g','LineWidth',1.4); grid on;
xlabel('Frequency (Hz)'); ylabel('Gain (dB)'); leg('aperture droop','equalizer 1/sinc(f\tau)');
title(sprintf('Aperture effect inside the band: %.4f dB droop at 50 Hz',droop50dB));
if saveFigs, savefig_('06_item3_recovery'); end

%% ---------------------------------------------------------------------
%  NUMBERS FOR THE SUMMARY TABLE
%  ---------------------------------------------------------------------
rms_ = @(x) sqrt(mean(x.^2));
fprintf('\n================ RESULTS (Digital Part A, items 1-3) ================\n');
fprintf('Nyquist rate = %d Hz ; 4x Nyquist = %d Hz\n',fNyq,4*fNyq);
fprintf('[Item 1] fs=100: max |recon-m| (ideal LPF, edge gain 1/2) = %.2e ; FIR check = %.2e\n',res1(1).maxErr,res1(1).maxErrFIR);
fprintf('[Item 1] fs=400: max |recon-m| (ideal LPF) = %.2e ; FIR check = %.2e\n',res1(2).maxErr,res1(2).maxErrFIR);
fprintf('[Item 1] fs=100, recovered 50 Hz amplitude for edge gain 1 / 0.5 / 0 = %.3f / %.3f / %.3f (true 0.5)\n',A50);
fprintf('[Item 1] sin(2pi50t) sampled at 100 Hz: max |sample| = %.1e (all zero)\n',sineMaxSample);
fprintf('[Item 2] fs=80: components after LPF at (Hz): %s ; RMS error vs original = %.3f V\n',mat2str(peaks),errAliasVsOrig);
fprintf('[Item 2] max |recon - (cos20 + 0.5cos30)| = %.2e  (output is exactly the aliased signal)\n',errAlias);
fprintf('[Item 3] Ts = %.3f ms, tau = %.4f ms (duty %.0f%%)\n',Ts*1e3,tau*1e3,d*100);
fprintf('[Item 3] sinc(50*tau) = %.5f -> droop %.4f dB (%.3f %%) ; sinc(20*tau) = %.5f\n',droop50,droop50dB,(1-droop50)*100,droop20);
fprintf('[Item 3] measured flat/natural line ratio: 20 Hz = %.5f, 50 Hz = %.5f (theory %.5f, %.5f)\n',ratio20,ratio50,droop20,droop50);
fprintf('[Item 3] equalizer boost needed at 50 Hz = +%.4f dB\n',boost50dB);
fprintf('[Item 3] for comparison, 100%% duty (tau=Ts): sinc(50*Ts) = %.4f -> %.3f dB\n',droopFull,droopFulldB);
fprintf('[Item 3] RMS recovery error: natural = %.2e ; flat (no eq) = %.2e ; flat + eq = %.2e\n',rms_(mr_nat-m),rms_(mr_flat_d-m),rms_(mr_eq-m));
fprintf('=====================================================================\n');

