function test_gpu_final()
    fprintf('=== 最终 GPU 测试 ===\n\n');
    
    % 1. 检查 .ptx 文件
    cuda_dir = 'C:\Users\Admin\Desktop\MMTools-main111\GMMNLSE\cuda';
    ptx_files = dir(fullfile(cuda_dir, '*.ptx'));
    fprintf('1. CUDA kernels:\n');
    fprintf('   找到 %d 个 .ptx 文件\n', length(ptx_files));
    for i = 1:length(ptx_files)
        fprintf('   ✅ %s (%.1f KB)\n', ptx_files(i).name, ptx_files(i).bytes/1024);
    end
    
    % 2. GPU 状态
    fprintf('\n2. GPU 状态:\n');
    gpu = gpuDevice;
    fprintf('   名称: %s\n', gpu.Name);
    fprintf('   可用显存: %.2f GB\n', gpu.AvailableMemory/1e9);
    
    % 3. GPU 计算测试（模拟 GMMNLSE 的操作）
    fprintf('\n3. GPU 计算性能测试:\n');
    
    N = 2^12;  % 和你的 Nt 一样
    num_modes = 4;  % λ₁ˣ, λ₁ʸ, λ₂ˣ, λ₂ʸ
    
    % CPU 测试
    fprintf('   CPU 测试...\n');
    tic;
    A_cpu = complex(rand(N, num_modes), rand(N, num_modes));
    for i = 1:10
        B_cpu = fft(A_cpu);
        C_cpu = B_cpu .* conj(B_cpu);
        A_cpu = ifft(C_cpu);
    end
    t_cpu = toc;
    
    % GPU 测试
    fprintf('   GPU 测试...\n');
    tic;
    A_gpu = gpuArray(complex(rand(N, num_modes), rand(N, num_modes)));
    for i = 1:10
        B_gpu = fft(A_gpu);
        C_gpu = B_gpu .* conj(B_gpu);
        A_gpu = ifft(C_gpu);
    end
    wait(gpu);
    t_gpu = toc;
    
    fprintf('\n   CPU 时间: %.3f 秒\n', t_cpu);
    fprintf('   GPU 时间: %.3f 秒\n', t_gpu);
    fprintf('   加速比: %.1fx\n', t_cpu/t_gpu);
    
    if t_gpu < t_cpu / 2
        fprintf('   ✅ GPU 加速有效！\n');
    else
        fprintf('   ⚠️  GPU 加速不明显\n');
    end
    
    clear A_cpu B_cpu C_cpu A_gpu B_gpu C_gpu;
    reset(gpu);
    
    fprintf('\n=== 测试完成 ===\n');
    fprintf('现在可以运行主代码了！记得设置 sim.cuda_dir_path\n');
end