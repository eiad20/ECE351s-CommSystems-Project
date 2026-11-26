N_bits = 64; #More bits for spectrum
bit_stream = randi([0,1], 1, N_bits)

dt = 0.1;
time_range = dt:dt:N_bits;

#Line Codes
AMI = [];
high = 1;

P_NRZ = [];

for (d = bit_stream)
  AMI = [AMI,  d*high*ones(1, 1/(2*dt)), zeros(1, 1/(2*dt))];
  high *= (-1)^d;

  P_NRZ = [P_NRZ, (2*d-1)*ones(1, 1/dt)];
endfor

figure('Position', [100, 100, 960, 360]);
plot(time_range, AMI)
hold on
plot(time_range, P_NRZ)
ylim([-1.2,1.2])
print("digB_1_figures/time_domain_waveform.png")
hold off;

df = 1/N_bits;
Fmax = 1/(2*dt);
freq_range = -Fmax:df:Fmax-df;

AMI_Spectrum = fftshift(fft(AMI));
P_NRZ_Spectrum = fftshift(fft(P_NRZ));

AMI_Power = (dt/length(freq_range)) * abs(AMI_Spectrum).^2;
P_NRZ_Power = (dt/length(freq_range)) * abs(P_NRZ_Spectrum).^2;

##figure('Position', [100, 100, 960, 360]);
##subplot(1,2,1);
##plot(freq_range, AMI_Power);
##subplot(1,2,2);
##plot(freq_range, P_NRZ_Power);
##print("digB_1_figures/Power_4096bits.png");
