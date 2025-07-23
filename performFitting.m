%PERFORMFITTING Multi-stage diode fitting
%   [OPT_PARAMS, FIT_RESULTS] = PERFORMFITTING(V, JD, PARAMS, CONFIG)
%   V and JD are vectors of measured voltage and current density.
%   PARAMS contains initial guesses, bounds and scaling factors.
%   CONFIG holds physical constants and optimisation settings.
%   The routine runs three stages:
%     1) fit_negative_region - fits Rsh and non-ohmic factor using negative V.
%     2) fit_positive_region - refines J0 and Rs in positive V ranges.
%     3) final_optimization  - global fit using LM/TR algorithms.
%   It returns the optimised parameter vector and modelled currents.
function [optimized_params, fit_results] = performFitting(data_V, data_JD, params, config)
    try
        % 缩放初始参数
        x0_scaled = params.x0 ./ params.scaleFactors;
        
        % 打印初始参数
        fprintf('\n初始参数（使用Lambert W函数估计）：\n');
        fprintf('J01 = %.6e A\n', params.x0(1));
        fprintf('Rs  = %.6e Ohm\n', params.x0(2));
        fprintf('Rsh = %.6e Ohm\n', params.x0(3));
        fprintf('k   = %.6e\n', params.x0(4));
        fprintf('J02 = %.6e A\n', params.x0(5));
        
        % 确保参数物理合理性
        if params.x0(2) <= 0  % Rs必须为正
            fprintf('警告: 初始Rs为负值或零，已自动调整为正值\n');
            params.x0(2) = max(params.lb(2), abs(params.x0(2)));
            x0_scaled = params.x0 ./ params.scaleFactors;
        end
        
        % 设置优化选项 - 尝试不同的算法
        fprintf('\n开始使用levenberg-marquardt算法进行拟合...\n');
        options_lm = optimoptions('lsqnonlin', ...
            'Display', 'iter-detailed', ...
            'Algorithm', 'levenberg-marquardt', ...
            'FunctionTolerance', 1e-8, ...
            'OptimalityTolerance', 1e-8, ...
            'StepTolerance', 1e-8, ...
            'MaxFunctionEvaluations', 8000, ...
            'FiniteDifferenceType', 'central', ... % 使用中心差分更精确
            'FiniteDifferenceStepSize', 1e-6, ... % 设置更合适的差分步长
            'DiffMaxChange', 1e-1, ...
            'DiffMinChange', 1e-8, ...
            'MaxIterations', 4000);
            
         % Combine optimisation options with parallel flag
        optcfg = config.optimization;
        if isfield(config, 'parallel') && isfield(config.parallel, 'use')
            optcfg.use_parallel = config.parallel.use;
        else
            optcfg.use_parallel = false;
        end

        % Initialize attempt counter
        attempt = 1;
        
       % 分阶段优化
        [x0_scaled, rel_errors] = fit_negative_region(data_V, data_JD, x0_scaled, params, config, options_lm, optcfg);
        [x0_scaled, rel_errors] = fit_positive_region(data_V, data_JD, x0_scaled, params, config, options_lm, rel_errors, optcfg);
        [optimized_params, fit_results] = final_optimization(data_V, data_JD, x0_scaled, params, config, options_lm, optcfg, attempt);
        
        % Check errors and adjust m in a small range when they are large
        rel = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
        if mean(rel) >= 2 || max(rel) >= 5
            fprintf('未达到误差阈值，尝试在%.1f到%.1f范围内调整m...\n', ...
                config.optimization.m_range(1), config.optimization.m_range(2));
            best_err = mean(rel);
            best_params = optimized_params;
            best_fit = fit_results;
            original_m = config.physics.m;
            m_values = config.optimization.m_range(1):config.optimization.m_step:config.optimization.m_range(2);
            for m_val = m_values
                config.physics.m = m_val;
                [tmp_params, tmp_fit] = final_optimization(data_V, data_JD, best_params ./ params.scaleFactors, params, config, options_lm, optcfg);
                tmp_rel = abs((tmp_fit.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
                tmp_err = mean(tmp_rel);
                if tmp_err < best_err
                    best_err = tmp_err;
                    best_params = tmp_params;
                    best_fit = tmp_fit;
                end
            end
            config.physics.m = original_m;
            optimized_params = best_params;
            fit_results = best_fit;
        end
    catch ME
        error('拟合过程出错: %s', ME.message);
    end
end

% 部分参数优化的误差函数
function err = errorFunctionNegative(x, data_V, data_JD, params, config, prior)
    % 反缩放参数
    x_actual = x .* params.scaleFactors;
    
    % 计算模型预测值
    predicted = diodeModel(data_V, x_actual, config);
    
    % 计算误差
    err = zeros(size(data_JD));
    
    for i = 1:length(data_JD)
        actual_abs = abs(data_JD(i));
        pred_abs = abs(predicted(i));
        
        threshold = 1e-12;
        
        if actual_abs < threshold || pred_abs < threshold
            err(i) = (predicted(i) - data_JD(i)) / max(1e-12, max(max(abs(data_JD))));
        else
            % 使用对数误差
            log_actual = log10(actual_abs);
            log_pred = log10(pred_abs);
            err(i) = log_pred - log_actual;
            
            % 保持符号一致性
            if sign(predicted(i)) ~= sign(data_JD(i))
                err(i) = err(i) * 4;
            end
        end
        
        % 针对更负的电压区域增加权重
        if data_V(i) < -0.3
            err(i) = err(i) * 3;
        end
    end

    if nargin < 6 || isempty(prior)
        prior = config.regularization.prior;
    end
    if isfield(config, 'regularization') && config.regularization.lambda > 0
        penalty = sqrt(config.regularization.lambda) * (x_actual(:) - prior(:));
        err = [err; penalty];
    end
end
