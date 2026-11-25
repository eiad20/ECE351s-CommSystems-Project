N_bits = 64;
bit_stream = randi([0,1], 1, N_bits)
#bit_stream = [1 0 1 0 1 0 1 0]
#bit_stream = ones(1, N_bits)
#bit_stream = zeros(1, N_bits)

dt = 0.1;
time_range = dt:dt:N_bits;

#Line Codes
AMI = [];
high = 1;

P_NRZ = [];

for (i = bit_stream)
  if (i==1)
    AMI = [AMI, high*ones(1, 1/(2*dt)), zeros(1, 1/(2*dt))];
    high *= -1;

    P_NRZ = [P_NRZ, ones(1, 1/dt)];
  else
    AMI = [AMI, zeros(1, 1/dt)];

    P_NRZ = [P_NRZ, -1*ones(1, 1/dt)];
  endif
endfor

##plot(time_range, AMI)
##hold on
##plot(time_range, P_NRZ)
##ylim([-1.2,1.2])
##print("digB_1_figures/time_domain_waveform.png")
##hold off;

df = 1/N_bits;
Fs = 1/(2*dt);
freq_range = -Fs:df:Fs-df;

AMI_Spectrum = fftshift(fft(AMI));
P_NRZ_Spectrum = fftshift(fft(P_NRZ));

##figure('Position', [100, 100, 960, 360]);
##subplot(1,2,1);
##plot(freq_range, AMI_Spectrum);
##subplot(1,2,2);
##plot(freq_range, P_NRZ_Spectrum);
##print("spectrum_1024bits.png")
