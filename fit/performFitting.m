function [optimized_params, fit_results] = performFitting(data_V, data_JD, params, config)
    try
        % 缩放初始参数
        x0_scaled = params.x0 ./ params.scaleFactors;
        
        guiLog('');
        guiLog('初始参数（使用Lambert W函数估计）：');
        guiLog(sprintf('J0 = %.6e A', params.x0(1)));
        guiLog(sprintf('Rs = %.6e Ohm', params.x0(2)));
        guiLog(sprintf('Rsh = %.6e Ohm', params.x0(3)));
        guiLog(sprintf('k = %.6e', params.x0(4)));
        
        % 确保参数物理合理性
        if params.x0(2) <= 0  % Rs必须为正
            guiLog('警告: 初始Rs为负值或零，已自动调整为正值');
            params.x0(2) = max(params.lb(2), abs(params.x0(2)));
            x0_scaled = params.x0 ./ params.scaleFactors;
        end
        
        % 设置优化选项 - 尝试不同的算法
        guiLog('');
        guiLog('开始使用levenberg-marquardt算法进行拟合...');
        options_lm = optimoptions('lsqnonlin', ...
            'Display', 'iter-detailed', ...
            'Algorithm', 'levenberg-marquardt', ...
            'FunctionTolerance', 1e-10, ...
            'OptimalityTolerance', 1e-10, ...
            'StepTolerance', 1e-10, ...
            'MaxFunctionEvaluations', 8000, ...
            'FiniteDifferenceType', 'central', ... % 使用中心差分更精确
            'FiniteDifferenceStepSize', 1e-6, ... % 设置更合适的差分步长
            'DiffMaxChange', 1e-1, ...
            'DiffMinChange', 1e-8, ...
            'MaxIterations', 4000);
            
       % 分阶段优化
        [x0_scaled, rel_errors] = fit_negative_region(data_V, data_JD, x0_scaled, params, config, options_lm);
        [x0_scaled, rel_errors] = fit_positive_region(data_V, data_JD, x0_scaled, params, config, options_lm, rel_errors);
        [optimized_params, fit_results] = final_optimization(data_V, data_JD, x0_scaled, params, config, options_lm, config.parallel);
    catch ME
        guiLog(sprintf('拟合过程出错: %s', ME.message));
        rethrow(ME);
    end
end
