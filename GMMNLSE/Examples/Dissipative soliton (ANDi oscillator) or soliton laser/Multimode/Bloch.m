% 假设仿真数据包含四个复光场：A1_x, A1_y, A2_x, A2_y
fields = struct();
fields.A1_x = complex(rand(1, Nt)); % 1528 nm 波长的 x 偏振态
fields.A1_y = complex(rand(1, Nt)); % 1528 nm 波长的 y 偏振态
fields.A2_x = complex(rand(1, Nt)); % 1559 nm 波长的 x 偏振态
fields.A2_y = complex(rand(1, Nt)); % 1559 nm 波长的 y 偏振态
% 定义窗函数（例如，汉宁窗）
window = hann(Nt);

% 对信号进行频域滤波
signal_fft = fft(fields.A1_x);  % FFT 变换
signal_fft_filtered = signal_fft .* window;  % 频域滤波
fields.A1_x_filtered = ifft(signal_fft_filtered);  % 反变换得到滤波后的信号
% 计算Stokes参数
I = abs(fields.A1_x).^2 + abs(fields.A1_y).^2 + abs(fields.A2_x).^2 + abs(fields.A2_y).^2; % 总光强
Q = abs(fields.A1_x).^2 + abs(fields.A2_x).^2 - abs(fields.A1_y).^2 - abs(fields.A2_y).^2; % 偏振光强度差
U = 2 * (real(conj(fields.A1_x) .* fields.A1_y + conj(fields.A2_x) .* fields.A2_y)); % 交叉项
V = 2 * (imag(conj(fields.A1_x) .* fields.A1_y - conj(fields.A2_x) .* fields.A2_y)); % 旋光项

% Stokes参数 S = [I, Q, U, V]
S = [I; Q; U; V]; % 每行代表一个 Stokes 参数
figure;
subplot(3, 1, 1);
plot(t, S(1, :)); % 绘制总光强 I
title('Stokes Parameter I');
xlabel('Time (ps)');
ylabel('Intensity (a.u.)');

subplot(3, 1, 2);
plot(t, S(2, :)); % 绘制 Q 参数
title('Stokes Parameter Q');
xlabel('Time (ps)');
ylabel('Polarization Difference (a.u.)');

subplot(3, 1, 3);
plot(t, S(3, :)); % 绘制 U 参数
title('Stokes Parameter U');
xlabel('Time (ps)');
ylabel('Cross Polarization (a.u.)');

