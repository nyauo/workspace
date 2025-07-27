function main()
    scriptDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(scriptDir, 'io')));
    addpath(genpath(fullfile(scriptDir, 'result')));
    addpath(genpath(fullfile(scriptDir, 'fit')));
    
    outputDir = 'outputs';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    % 加载配置和数据
    config = loadConfig();
    [data_V, data_JD] = loadData();
    
    % 根据正向数据估计理想因子并覆盖配置
    n_est = estimateIdealityFactor(data_V, data_JD, config);
    fprintf('估计的理想因子 n = %.3f\n', n_est);
    config.physics.n = n_est;
    
    % 验证输入数据
    validateInputData(data_V, data_JD);
    
    % 询问是否从历史参数文件加载初始参数
    use_history = input('是否使用历史参数文件作为初始参数? (y/n): ', 's');
    if strcmpi(use_history, 'y')
        % 列出可用的参数文件
        mat_files = dir(fullfile(outputDir, 'adjusted_params_*.mat'));
        txt_files = dir(fullfile(outputDir, 'adjusted_params_*.txt'));
        
        if isempty(mat_files)
            fprintf('未找到历史参数文件，将使用Lambert W函数估计初始参数\n');
            params = initializeParameters(data_V, data_JD, config);
        else
            % 显示可用文件
            fprintf('可用的参数文件:\n');
            for i = 1:length(mat_files)
                fprintf('%d: %s\n', i, mat_files(i).name);
            end
            
            % 让用户选择文件
            file_idx = input('请选择要加载的文件编号(输入0取消): ');
            if file_idx > 0 && file_idx <= length(mat_files)
                % 加载选定的文件
                load_file = fullfile(outputDir, mat_files(file_idx).name);
                loaded_data = load(load_file);
                
                % 提取参数
                if isfield(loaded_data, 'params')
                    fprintf('从文件 %s 加载参数\n', load_file);
                    
                    % 创建参数结构体
                    params = struct();
                    params.x0 = loaded_data.params;
                    
                    % 参数范围设置
                    params.ub = [1e-6, 1e4, 1e10, 1e-5];    % 上界
                    params.lb = [1e-12, 1e1, 1e5, 1e-10];    % 下界
                    
                    % 确保初始值在范围内
                    params.x0 = min(max(params.x0, params.lb), params.ub);
                    
                    % 缩放因子
                    params.scaleFactors = [1e-9, 1e3, 1e7, 1e-8];
                    
                    % 显示加载的参数
                    fprintf('加载的参数:\n');
                    fprintf('J0 = %.6e A\n', params.x0(1));
                    fprintf('Rs = %.6e Ohm\n', params.x0(2));
                    fprintf('Rsh = %.6e Ohm\n', params.x0(3));
                    fprintf('k = %.6e\n', params.x0(4));
                else
                    fprintf('文件格式错误，将使用Lambert W函数估计初始参数\n');
                    params = initializeParameters(data_V, data_JD, config);
                end
            else
                fprintf('取消加载参数文件，将使用Lambert W函数估计初始参数\n');
                params = initializeParameters(data_V, data_JD, config);
            end
        end
    else
        % 初始化参数 - 传入数据，用于Lambert W函数估计初始参数
        params = initializeParameters(data_V, data_JD, config);
    end
    
    % 执行拟合，performFitting 内部会按电压区域调用子函数逐步优化
    [optimized_params, fit_results] = performFitting(data_V, data_JD, params, config);
    
    % 计算各分量电流
    currents = calculateCurrents(data_V, optimized_params, config);
    
    % 绘制结果
    plotResults(data_V, data_JD, fit_results, currents);
    
    % 保存拟合结果和图像
    save_results = input('是否保存拟合结果? (y/n): ', 's');
    if strcmpi(save_results, 'y')
        saveResults(data_V, data_JD, optimized_params, fit_results, currents);
    end
    
    % 询问用户是否满意拟合结果，如果不满意则进入交互式调整模式
    interactive_adjust = input('是否进入交互式参数调整模式? (y/n): ', 's');
    if strcmpi(interactive_adjust, 'y')
        [refined_params, refined_fit] = interactiveParameterAdjustment(data_V, data_JD, optimized_params, config);
        
        % 计算新的分量电流
        refined_currents = calculateCurrents(data_V, refined_params, config);
        
        % 绘制新的结果
        figure;
        plotResults(data_V, data_JD, refined_fit, refined_currents);
        

        
        % 保存调整后的参数
        saveAdjustedParameters(refined_params);
    end
end

