% clear;
% close all;
% 
% % 假设你已经从仿真或测量中获得了 Stokes 向量的数据
% % 这里使用示例数据代替
% Nt = 1000;  % 时间点数
% t = linspace(0, 100, Nt);  % 假设时间范围
% 
% % 示例 Stokes 向量数据，I, Q, U, V
% I = rand(1, Nt) + 1;  % 总光强
% Q = sin(t);           % 偏振光强度差
% U = cos(t);           % 交叉偏振
% V = sin(2*t);         % 旋光项
% 
% % 计算归一化的 Sx, Sy, Sz
% Sx = Q ./ I;  % Q/I
% Sy = U ./ I;  % U/I
% Sz = V ./ I;  % V/I
% 
% norm_S = sqrt(Sx.^2 + Sy.^2 + Sz.^2);
% Sx = Sx ./ norm_S;
% Sy = Sy ./ norm_S;
% Sz = Sz ./ norm_S;
% 
% 
% Sx_smooth = smooth(Sx, 0.1, 'loess');
% Sy_smooth = smooth(Sy, 0.1, 'loess');
% Sz_smooth = smooth(Sz, 0.1, 'loess');
% 
% % 绘制 Bloch 球
% figure;
% 
% % 使用 meshgrid 创建单位球的坐标
% [phi, theta] = meshgrid(linspace(0, 2*pi, 100), linspace(0, pi, 50));
% x = sin(theta) .* cos(phi);
% y = sin(theta) .* sin(phi);
% z = cos(theta);
% 
% % 绘制球面
% surf(x, y, z, 'EdgeColor', 'none', 'FaceAlpha', 0.1);
% hold on;
% 
% % 绘制 Stokes 向量的轨迹
% plot3(Sx, Sy, Sz, 'LineWidth', 2, 'Color', 'b');
% 
% % 添加标签和标题
% title('Stokes Vector Trajectory on Bloch Sphere');
% xlabel('S_x');
% ylabel('S_y');
% zlabel('S_z');
% axis equal;
% grid on;
% 生成双色 Stokes 向量数据
t = linspace(0, 2*pi, 200); 
Sx1 = 0.1; 
Sy1 = cos(t); 
Sz1 = 0.5 * sin(2*t); 

Sx2 = 0.2; 
Sy2 = cos(t + pi/0.1); 
Sz2 = 0.5 * sin(2*t + pi/5);

% 归一化 Stokes 向量
norm_S1 = sqrt(Sx1.^2 + Sy1.^2 + Sz1.^2);
norm_S2 = sqrt(Sx2.^2 + Sy2.^2 + Sz2.^2);

Sx1 = Sx1 ./ norm_S1;
Sy1 = Sy1 ./ norm_S1;
Sz1 = Sz1 ./ norm_S1;

Sx2 = Sx2 ./ norm_S2;
Sy2 = Sy2 ./ norm_S2;
Sz2 = Sz2 ./ norm_S2;

% 绘制 Bloch 球
figure;

% 创建单位球
[phi, theta] = meshgrid(linspace(0, 2*pi, 100), linspace(0, pi, 50));
x = sin(theta) .* cos(phi);
y = sin(theta) .* sin(phi);
z = cos(theta);

% 绘制球面
surf(x, y, z, 'EdgeColor', 'none', 'FaceAlpha', 0.1);
hold on;

% 绘制 Stokes 向量轨迹（双色）
plot3(Sx1, Sy1, Sz1, 'LineWidth', 2, 'Color', 'b'); % 第一波长
plot3(Sx2, Sy2, Sz2, 'LineWidth', 2, 'Color', 'r'); % 第二波长

% 设置图形
title('双色 Stokes Vector Trajectory on Bloch Sphere');
xlabel('S_x');
ylabel('S_y');
zlabel('S_z');
axis equal;
grid on;
