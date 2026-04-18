% 能量管理孤子

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');
    
%% Gain info
gain_rate_eqn.gain_medium = 'Er'; % 指定增益介质
gain_rate_eqn.base_medium = 'silica'; % 指定基本介质
gain_rate_eqn.reuse_data = false; % 对于环形或线性空腔，脉冲最终会进入稳定状态。
                                  % 如果重复使用上一次往返的泵和 ASE 数据，收敛速度会更快，特别是在反回泵时。
gain_rate_eqn.linear_oscillator = false; % 对于线性振荡器来说，两个方向的脉冲会同时出现，从而耗尽增益；
                                         % 因此，需要考虑后向传播脉冲。
gain_rate_eqn.core_diameter = 8; % um 纤芯直径
gain_rate_eqn.cladding_diameter = gain_rate_eqn.core_diameter; % um 包层直径
gain_rate_eqn.core_NA = 0.13;  % 数值孔径
gain_rate_eqn.absorption_wavelength_to_get_N_total = 1530; % nm
gain_rate_eqn.absorption_to_get_N_total = 80; % dB/m
gain_rate_eqn.pump_wavelength = 1480; % nm 泵浦波长
gain_rate_eqn.copump_power = 0.100; % W 泵浦功率
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % 在 LMA 光纤中，由于信号场的 ASE 模式数可能大于 1，因此使用该系数来正确考虑 ASE。 如果是 [] 这样的空值，则为 length(sim.midx)。
gain_rate_eqn.max_iterations = 0; % 对于反回泵或考虑 ASE，需要进行迭代。
gain_rate_eqn.tol = 1e-5; % 迭代的宽容度
gain_rate_eqn.verbose = true; % 在计算增益的迭代过程中显示信息（最终脉冲能量）

%% 设置光纤参数
% 一般参数
sim.lambda0 = 1550e-9; % m; 中心波长
sim.f0 = 2.99792458e-4/sim.lambda0; % THz; 中心频率
%sim.progress_bar = false;
sim.gpu_yes = false; % 是否使用 GPU

sim_ND = sim; % 无源光纤
sim_ND.progress_bar_name = 'SMF (10.4um)';
sim_ND.save_period = 0.1; % m 每隔多少距离保存一次完整的光场数据

[fiber_ND,sim_ND] = load_default_GMMNLSE_propagate([],sim_ND,'single_mode'); % 无源光纤
fiber_ND.MFD = 10.4;

% -------------------------------------------------------------------------
% Normal-dpsersion SMF
fiber_ND.L0 = 0.9;
fiber_ND.betas = [5.8126e6, 4.89e3, -22.8e-3]';% ps^k/m
fiber_ND1 = fiber_ND;
fiber_ND1.L0 = 2.6;
% fiber_ND1.betas = [5.8126e6, 4.89e3, -22.8e-3, 0.11e-3, -8.5e-7]';% ps^k/m

% 增益光纤
sim_Gain = sim;
sim_Gain.gain_model = 2; % 使用速率-增益模型
sim_Gain.progress_bar_name = 'Gain (9.5um)';
sim_Gain.save_period = 0.01; % m 每隔多少距离保存一次完整的光场数据
fiber_Gain.L0 = 1.2;
fiber_Gain.MFD = 9.5;
fiber_Gain.betas = [5.8126e6, 5.0e3, -25.5e-3]';% ps^k/m
% fiber_Gain.betas = [5.8126e6, 5.0e3, -25.5e-3, 0.13e-3, -6.5e-7]';% ps^k/m
[fiber_Gain,sim_Gain] = load_default_GMMNLSE_propagate(fiber_Gain,sim_Gain,'single_mode'); % 对于增益光纤

SA_length = 0.5; % m
SF_length = 0.5; % m
component_z = [SA_length, SF_length];

% -------------------------------------------------------------------------
% ----------------------------------- All ---------------------------------
% -------------------------------------------------------------------------
fiber_cavity = [fiber_ND fiber_Gain fiber_ND1];
sim_cavity = [sim_ND sim_Gain sim_ND];

%% 设置一般空腔参数
max_rt = 500; % 最大往返次数（以防不收敛）
Nt = 2^12; % 点数
time_window = 100; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm
OC = 0.9; % 输出耦合
loss = 0; % 腔体总损耗
saturation_power = 100; % 可饱和吸收体的饱和功率；W
moddepth = 0.35; % 可饱和吸收体的调制深度
tol_convergence = 1e-3; % 输出脉冲能量收敛的容差

%% 频谱滤波器参数
gaussexpo = 1;
plot_filter_yes = true;
spectral_filter = struct('bw',10, ... % 带宽 (nm)
                         'cw',1569); % 中心波长 (nm)

%% 设置初始条件
tfwhm = 0.5; % ps
total_energy = 1; % nJ

prop_output = build_MMgaussian(tfwhm, time_window, total_energy,1,Nt);

%% 已保存的光场信息
L0 = sum([fiber_cavity.L0]); % 光纤总长度
save_num = sum(int64([fiber_cavity.L0]./[sim_cavity.save_period])) + 1; % 要保存的光场总数
save_num = double(save_num);
saved_z = zeros(1,save_num); % 传播距离
field = cell(1,max_rt); 
splice_z = cumsum([fiber_ND.L0, fiber_Gain.L0, SA_length, SF_length, fiber_ND1.L0]);% 接合点的位置
output_field = zeros(Nt,1,max_rt); % 输出脉冲
max_save_per_fiber = 30;

%% 加载增益参数
L_air = 0.5; % 1 是自由空间长度
c = 299792458; % m/s
v = 1/fiber_cavity(1).betas(2)*1e12; % 在光纤中传播的速度

t_rep = L0/v + L_air/c; % s; 完成一次往返所需的时间（脉冲的反重复率）
                        % 该增益模型求解的是稳态条件下光纤的增益；因此，与掺杂离子的寿命相比，重复率必须很高。
gain_rate_eqn.t_rep = t_rep;

gain_rate_eqn = gain_info( fiber_Gain,sim_Gain,gain_rate_eqn,ifftshift(lambda,1) );

% %% 通过空间滤波器的传输矩阵
% spatial_hwhm = 5e-6; % um
% filter_type = 'Gaussian';
% offcenter = [0,0]*1e-6; % m
% spatial_filter_matrix = calc_filter_matrix(1, 1*1e-6, spatial_hwhm, offcenter, filter_type);

%% 运行空腔模拟

func = analyze_sim; % 若干分析功能的容器
% 初始化一些参数
output_energy = zeros(max_rt,1);
rt_num = 0;
pulse_survives = true;
while rt_num < max_rt
    time_delay = 0;
    current_z = 0;
    zn = 1;
    rt_num = rt_num + 1;
    
    t_iteration_start = tic;
    cprintf('*[1 0.5 0.31]','Iteration %d', rt_num);
    % -----------------------------------------------------------------
    % SMF1
    prop_output = GMMNLSE_propagate(fiber_cavity(1), prop_output, sim_cavity(1));
    time_delay = time_delay + prop_output.t_delay(end);
    
    [saved_field,saved_z_this_fiber] = func.extract_saved_field(prop_output.fields,max_save_per_fiber,current_z,prop_output.z);
    field{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_field;
    saved_z(zn:zn+size(saved_field,3)-1) = saved_z_this_fiber;

    current_z = saved_z_this_fiber(end);
    zn = zn + size(saved_field,3)-1;
    
    % -----------------------------------------------------------------
    % EDF
    prop_output = GMMNLSE_propagate(fiber_cavity(2), prop_output, sim_cavity(2), gain_rate_eqn);
    time_delay = time_delay + prop_output.t_delay(end);
   
    [saved_field,saved_z_this_fiber] = func.extract_saved_field(prop_output.fields,max_save_per_fiber,current_z,prop_output.z);
    field{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_field;
    saved_z(zn:zn+size(saved_field,3)-1) = saved_z_this_fiber;

    current_z = saved_z_this_fiber(end);
    zn = zn + size(saved_field,3);

    % % -----------------------------------------------------------------
    % % 空间滤波器
    % prop_output = spatial_filter_moderesolved(prop_output,spatial_filter_matrix);
    
    % -----------------------------------------------------------------
    % 输出耦合器
    output_field(:,:,rt_num) = sqrt(OC)*prop_output.fields(:,:,end);
    prop_output.fields = sqrt(1-OC)*prop_output.fields(:,:,end);

    % -----------------------------------------------------------------
    % 饱和吸收体
    prop_output = saturable_absorber_action_simple(prop_output, saturation_power, moddepth);
    [saved_field,saved_z_SA] = func.extract_saved_field(prop_output.fields, 1, current_z, SA_length);

    field{rt_num}(:,:,zn) = saved_field;
    saved_z(zn) = saved_z_SA; 

    current_z = current_z + SA_length;
    zn = zn + size(saved_field, 3);
    % -----------------------------------------------------------------
    % 若输出能量，则把下面三行注释掉
    field{rt_num}(:,:,zn) = saved_field;
    saved_z(zn) = current_z;
    zn = zn + 1;
    % -----------------------------------------------------------------
    % 频谱滤波器
    if rt_num ~= 1
        close(fig_filter); % 关闭上一图，绘制新图
    end
    [prop_output,fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter.cw, spectral_filter.bw, gaussexpo ,plot_filter_yes); % Filter
    prop_output.fields = sqrt(1-loss)*prop_output.fields(:,:,end);
    [saved_field,saved_z_SF] = func.extract_saved_field(prop_output.fields, 1, current_z, SF_length);

    field{rt_num}(:,:,zn) = saved_field;
    saved_z(zn) = saved_z_SF; % 末端位置

    current_z = current_z + SF_length;
    zn = zn + size(saved_field, 3);

    % -----------------------------------------------------------------
    % SMF2
    prop_output = GMMNLSE_propagate(fiber_cavity(3), prop_output, sim_cavity(3));
    time_delay = time_delay + prop_output.t_delay(end);
   
    [saved_field,saved_z_this_fiber] = func.extract_saved_field(prop_output.fields,max_save_per_fiber,current_z,prop_output.z);
    field{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_field;
    saved_z(zn:zn+size(saved_field,3)-1) = saved_z_this_fiber;

    current_z = saved_z_this_fiber(end);
    zn = zn + size(saved_field,3)-1;

    saved_z = saved_z(1:zn);
    % -----------------------------------------------------------------
    % 输出场的能量
    output_energy(rt_num) = trapz(abs(output_field(:,:,rt_num)).^2)*prop_output.dt/1e3; % 能量（nj)

    % 如果能量停止变化，我们就停止了！
    if rt_num ~= 1
        close(fig); % 关闭上一图，绘制新图
    end
    warning('off')
    output_field(:,:,rt_num) = pulse_tracker(output_field(:,:,rt_num));
    warning('on');
    [converged_yes,fig] = check_convergence( output_energy,output_field(:,:,rt_num),f,t,tol_convergence,true );
    
    % ---------------------------------------------------------------------
    % 显示每次往返的运行时间
    t_iteration_end = toc(t_iteration_start);
    t_iteration_spent = datevec(t_iteration_end/3600/24);
    fprintf(': %1u:%2u:%3.1f\n',t_iteration_spent(4),t_iteration_spent(5),t_iteration_spent(6));
    
    % ---------------------------------------------------------------------
    % 根据 "time_delay "更新重复率
    % 脉冲相对于运动帧会发生时间上的偏移。
    % 因为我在代码中实现了脉冲居中功能，我可以使用每次往返中的 "time_delay "信息进行校准，从而得到实际的重复率。
    gain_rate_eqn.t_rep = t_rep + time_delay*1e-12;
    
    % ------------------------ FROG、沿光纤的光场信息 --------------------------------
    if rt_num ~= 1
        close(fig_evolution); % 关闭上一图，绘制新图
        close(fig_FROG);
    end
    fig_evolution = func.analyze_fields(t,f,field{rt_num},saved_z,splice_z);
    fig_FROG = func.analyze_FROG(t,f,output_field(:,:,rt_num), 'SHG-FROG');

    % ------------------------------- 自相关 ----------------------------------------
    if rt_num ~= 1
        close(fig_autocorr); % 关闭上一图，绘制新图
    end
    % 计算自相关
    G1 = func.autocorr1(output_field(:,:,rt_num));
    G2 = func.autocorr2(output_field(:,:,rt_num));

    % 可视化并诊断脉冲质量
    fig_autocorr = func.analyze_autocorr(t, G1, G2, rt_num, tfwhm);

    % ---------------------------------------------------------------------
    % Break if converged
    if converged_yes
        cprintf('blue','The field has converged!\n');
        break;
    end
    % Break if pulse dies
    if output_energy(rt_num) < 0.001 % 小于 0.001 nJ
        disp('The pulse dies.');
        pulse_survives = false;
        break;
    end
end

%% 完成模拟并保存数据
% 清除数据中的缩减部分
field = field(1:rt_num);
output_field = output_field(:,:,1:rt_num);
energy = output_energy(arrayfun(@any,output_energy)); 

% close(fig,fig_filter,fig_evolution);

%% Compress the pulse 压缩脉冲

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,output_field(:,:,end),'Treacy-t',pi/6,1e-3/600,true );