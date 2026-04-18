function test_gmmnlse_cuda()
    % 找到 GMMNLSE_propagate 的位置
    gmmnlse_file = which('GMMNLSE_propagate');
    if isempty(gmmnlse_file)
        error('找不到 GMMNLSE_propagate.m，请检查路径');
    end
    
    fprintf('GMMNLSE_propagate 位置:\n  %s\n\n', gmmnlse_file);
    
    % 推测 cuda 目录
    gmmnlse_dir = fileparts(gmmnlse_file);
    cuda_dir = fullfile(gmmnlse_dir, 'cuda');
    
    fprintf('预期的 cuda 目录:\n  %s\n\n', cuda_dir);
    
    % 检查目录是否存在
    if exist(cuda_dir, 'dir')
        fprintf('✅ cuda 目录存在\n');
        
        % 列出 .cu 文件
        cu_files = dir(fullfile(cuda_dir, '*.cu'));
        fprintf('找到 %d 个 .cu 文件:\n', length(cu_files));
        for i = 1:min(5, length(cu_files))
            fprintf('  %s\n', cu_files(i).name);
        end
        
        % 列出 .ptx 文件
        ptx_files = dir(fullfile(cuda_dir, '*.ptx'));
        fprintf('\n找到 %d 个 .ptx 文件:\n', length(ptx_files));
        for i = 1:min(5, length(ptx_files))
            fprintf('  %s\n', ptx_files(i).name);
        end
        
        if isempty(ptx_files)
            fprintf('\n❌ 没有 .ptx 文件，需要编译！\n');
            fprintf('运行以下命令手动编译:\n');
            fprintf('  cd(''%s'');\n', cuda_dir);
            fprintf('  !nvcc -ptx GMMNLSE_nonlinear_sum_with_polarization.cu\n');
        end
    else
        fprintf('❌ cuda 目录不存在！\n');
    end
end