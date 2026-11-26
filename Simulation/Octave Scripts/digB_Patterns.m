N_bits = 8;
bit_stream = [1 0 1 0 1 0 1 0]
bit_stream = [bit_stream; ones(1, N_bits)]
bit_stream = [bit_stream; zeros(1, N_bits)]

dt = 0.1;
time_range = dt:dt:N_bits;

#Line Codes
AMI = ones(3, N_bits);
high = 1;

P_NRZ = [3, N_bits];


for (n = 1:3)
  count = 1;
  for (d = bit_stream(n, 1:N_bits))
    AMI(n, count:(count + 1/dt - 1)) = [d*high*ones(1, 1/(2*dt)), zeros(1, 1/(2*dt))];
    high *= (-1)^d;

    P_NRZ(n, count:(count + 1/dt - 1)) = [(2*d-1)*ones(1, 1/dt)];
    count = count + 1/dt;
  endfor
endfor


figure('Position', [100, 100, 1440, 360]);
subplot(1,3,1);
plot(time_range, AMI(1, :))
hold on
plot(time_range, P_NRZ(1, :))
ylim([-1.2,1.2])

subplot(1,3,2);
plot(time_range, AMI(2, :))
hold on
plot(time_range, P_NRZ(2, :))
ylim([-1.2,1.2])

subplot(1,3,3);
plot(time_range, AMI(3, :))
hold on
plot(time_range, P_NRZ(3, :))
ylim([-1.2,1.2])
#print("digB_1_figures/CLK_Recovery_Patterns.png")
hold off;
