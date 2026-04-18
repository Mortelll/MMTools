function diagnose_gpu()
    fprintf('=== GPU 使用诊断 ===\n\n');
    
    % 1. 检查 GPU 状态
    fprintf('1. GPU 状态:\n');
    gpu = gpuDevice;
    fprintf('   名称: %s\n', gpu.Name);
    fprintf('   显存: %.2f GB 可用 / %.2f GB 总共\n', ...
        gpu.AvailableMemory/1e9, gpu.TotalMemory/1e9);
    
    % 2. 检查 CUDA 路径
    fprintf('\n2. CUDA 环境:\n');
    cuda_path = getenv('CUDA_PATH');
    if isempty(cuda_path)
        fprintf('   ❌ CUDA_PATH 未设置\n');
    else
        fprintf('   ✅ CUDA_PATH = %s\n', cuda_path);
    end
    
    % 3. 检查 .ptx 文件
    fprintf('\n3. CUDA kernel 编译状态:\n');
    cuda_dir = '../../../GMMNLSE algorithm/cuda/';
    ptx_files = dir([cuda_dir '*.ptx']);
    if isempty(ptx_files)
        fprintf('   ❌ 没有找到 .ptx 文件（未编译）\n');
    else
        fprintf('   ✅ 找到 %d 个 .ptx 文件\n', length(ptx_files));
        for i = 1:length(ptx_files)
            fprintf('      - %s\n', ptx_files(i).name);
        end
    end
    
    % 4. 测试 GPU 计算
    fprintf('\n4. GPU 计算测试:\n');
    N = 2^12;
    
    % CPU 测试
    tic;
    A_cpu = rand(N) + 1i*rand(N);
    B_cpu = fft(A_cpu);
    t_cpu = toc;
    
    % GPU 测试
    tic;
    A_gpu = gpuArray(rand(N) + 1i*rand(N));
    B_gpu = fft(A_gpu);
    wait(gpu);
    t_gpu = toc;
    
    fprintf('   CPU 时间: %.4f 秒\n', t_cpu);
    fprintf('   GPU 时间: %.4f 秒\n', t_gpu);
    fprintf('   加速比: %.1fx\n', t_cpu/t_gpu);
    
    if t_gpu < t_cpu
        fprintf('   ✅ GPU 加速有效\n');
    else
        fprintf('   ❌ GPU 没有加速（可能在用 CPU）\n');
    end
    
    clear A_cpu B_cpu A_gpu B_gpu;
    reset(gpu);
    
    fprintf('\n=== 诊断完成 ===\n');
end