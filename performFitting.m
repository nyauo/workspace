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

        % Initialize retry variables
        attempt = 1;
        retry_count = 0;
        
       % 分阶段优化
        [x0_scaled, rel_errors] = fit_negative_region(data_V, data_JD, x0_scaled, params, config, options_lm, optcfg);
        [x0_scaled, rel_errors] = fit_positive_region(data_V, data_JD, x0_scaled, params, config, options_lm, rel_errors, optcfg);
        [optimized_params, fit_results] = final_optimization(data_V, data_JD, x0_scaled, params, config, options_lm, optcfg, attempt, retry_count);
        
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

function [x0_scaled, rel_errors] = fit_negative_region(data_V, data_JD, x0_scaled, params, config, options_lm, optcfg)
    fprintf('\n第一阶段：优化Rsh和非欧姆系数k...\n');
    neg_idx = find(data_V < -0.1);
    if ~isempty(neg_idx)
        neg_V = data_V(neg_idx);
        neg_JD = data_JD(neg_idx);
        x0_limited = x0_scaled;
        param_mask = [false, false, true, true, false];
        neg_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_limited, param_mask, neg_V, neg_JD, params, config, config.regularization.prior, []);
        x0_neg_opt = x0_scaled(param_mask);
        lb_neg = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_neg = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        optcfg_neg = optcfg;
        [x_neg_opt, ~] = runWithMultiStart(neg_errFun, x0_neg_opt, lb_neg, ub_neg, options_lm, optcfg_neg);
        x0_scaled(param_mask) = x_neg_opt;
        x_actual_neg = x0_scaled .* params.scaleFactors;
        fit_JD_neg = diodeModel(neg_V, x_actual_neg, config);
        neg_rel_errors = abs((fit_JD_neg - neg_JD) ./ (abs(neg_JD) + eps)) * 100;
        fprintf('负电压区域拟合：平均相对误差 = %.2f%%, 最大相对误差 = %.2f%%\n', ...
            mean(neg_rel_errors), max(neg_rel_errors));
        fprintf('优化后的Rsh = %.6e Ohm\n', x_actual_neg(3));
        fprintf('优化后的k = %.6e\n', x_actual_neg(4));
    end
    % 计算整个数据集的相对误差供下一阶段使用
    x_all = x0_scaled .* params.scaleFactors;
    fit_all = diodeModel(data_V, x_all, config);
    rel_errors = abs((fit_all - data_JD) ./ (abs(data_JD) + eps));
end

function [x0_scaled, rel_errors] = fit_positive_region(data_V, data_JD, x0_scaled, params, config, options_lm, prev_errors, optcfg)
    fprintf('\n第二阶段：细分正电压区域优化...\n');
    low_pos_idx = find(data_V > 0 & data_V <= 0.15);
    high_pos_idx = find(data_V > 0.15);

    if ~isempty(low_pos_idx)
        fprintf('优化低正电压区域 (0 - 0.15V)...\n');
        low_pos_V = data_V(low_pos_idx);
        low_pos_JD = data_JD(low_pos_idx);
        param_mask = [true, false, false, false, false];
        low_pos_prev = prev_errors(low_pos_idx);
        low_pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, low_pos_V, low_pos_JD, params, config, config.regularization.prior, low_pos_prev);
        x0_low_pos_opt = x0_scaled(param_mask);
        lb_low_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_low_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        optcfg_low = optcfg;
        [x_low_pos_opt, ~] = runWithMultiStart(low_pos_errFun, x0_low_pos_opt, lb_low_pos, ub_low_pos, options_lm, optcfg_low);
        x0_scaled(param_mask) = x_low_pos_opt;
        % 更新误差供下一小阶段使用
        x_all = x0_scaled .* params.scaleFactors;
        fit_all = diodeModel(data_V, x_all, config);
        prev_errors = abs((fit_all - data_JD) ./ (abs(data_JD) + eps));
        x_actual_low_pos = x0_scaled .* params.scaleFactors;
        fit_JD_low_pos = diodeModel(low_pos_V, x_actual_low_pos, config);
        low_pos_rel_errors = abs((fit_JD_low_pos - low_pos_JD) ./ (abs(low_pos_JD) + eps)) * 100;
        fprintf('低正电压区域拟合：平均相对误差 = %.2f%%, 最大相对误差 = %.2f%%\n', ...
            mean(low_pos_rel_errors), max(low_pos_rel_errors));
        fprintf('优化后的J0 = %.6e A\n', x_actual_low_pos(1));
    end

    if ~isempty(high_pos_idx)
        fprintf('优化高正电压区域 (>0.15V)...\n');
        high_pos_V = data_V(high_pos_idx);
        high_pos_JD = data_JD(high_pos_idx);
        param_mask = [false, true, false, false, false];
        high_pos_prev = prev_errors(high_pos_idx);
        high_pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, high_pos_V, high_pos_JD, params, config, config.regularization.prior, high_pos_prev);
        x0_high_pos_opt = x0_scaled(param_mask);
        lb_high_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_high_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        optcfg_high = optcfg;
        [x_high_pos_opt, ~] = runWithMultiStart(high_pos_errFun, x0_high_pos_opt, lb_high_pos, ub_high_pos, options_lm, optcfg_high);
        x0_scaled(param_mask) = x_high_pos_opt;
        % 更新误差供下一小阶段使用
        x_all = x0_scaled .* params.scaleFactors;
        fit_all = diodeModel(data_V, x_all, config);
        prev_errors = abs((fit_all - data_JD) ./ (abs(data_JD) + eps));
        
        
        x_actual_high_pos = x0_scaled .* params.scaleFactors;
        fit_JD_high_pos = diodeModel(high_pos_V, x_actual_high_pos, config);
        high_pos_rel_errors = abs((fit_JD_high_pos - high_pos_JD) ./ (abs(high_pos_JD) + eps)) * 100;
        fprintf('高正电压区域拟合：平均相对误差 = %.2f%%, 最大相对误差 = %.2f%%\n', ...
            mean(high_pos_rel_errors), max(high_pos_rel_errors));
        fprintf('优化后的Rs = %.6e Ohm\n', x_actual_high_pos(2));
    end

    pos_idx = find(data_V > 0);
    if ~isempty(pos_idx)
        fprintf('\n综合优化正电压区域...\n');
        pos_V = data_V(pos_idx);
        pos_JD = data_JD(pos_idx);
        param_mask = [true, true, false, false, false];
        pos_prev = prev_errors(pos_idx);
        pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, pos_V, pos_JD, params, config, config.regularization.prior, pos_prev);
        x0_pos_opt = x0_scaled(param_mask);
        lb_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        optcfg_pos = optcfg;
        [x_pos_opt, ~] = runWithMultiStart(pos_errFun, x0_pos_opt, lb_pos, ub_pos, options_lm, optcfg_pos);
        x0_scaled(param_mask) = x_pos_opt;

        x_actual_pos = x0_scaled .* params.scaleFactors;
        fit_JD_pos = diodeModel(pos_V, x_actual_pos, config);
        pos_rel_errors = abs((fit_JD_pos - pos_JD) ./ (abs(pos_JD) + eps)) * 100;
        fprintf('正电压区域综合拟合：平均相对误差 = %.2f%%, 最大相对误差 = %.2f%%\n', ...
            mean(pos_rel_errors), max(pos_rel_errors));
        fprintf('优化后的J0 = %.6e A\n', x_actual_pos(1));
        fprintf('优化后的Rs = %.6e Ohm\n', x_actual_pos(2));
    end
    % 计算整个数据集的相对误差供下一阶段使用
    x_all = x0_scaled .* params.scaleFactors;
    fit_all = diodeModel(data_V, x_all, config);
    rel_errors = abs((fit_all - data_JD) ./ (abs(data_JD) + eps));    
end

function [optimized_params, fit_results] = final_optimization(data_V, data_JD, x0_scaled, params, config, options_lm, optcfg, attempt, retry_count)
    if nargin < 8
        attempt = 1;
    end
    if nargin < 9
        retry_count = 0;
    end
    fprintf('\n第三阶段：全区域拟合...\n');
    errFun = @(x) errorFunction(x, data_V, data_JD, params, config, config.regularization.prior);
    [x_scaled_optimized,resnorm_lm]=runWithMultiStart(errFun,x0_scaled,params.lb./ params.scaleFactors,params.ub ./ params.scaleFactors,options_lm,optcfg);
    residual_lm = []; exitflag_lm = []; output_lm = [];
    if x_scaled_optimized(2) * params.scaleFactors(2) <= 0
        fprintf('警告: LM算法产生了负值或零的Rs，正在调整为正值\n');
        x_scaled_optimized(2) = params.lb(2) / params.scaleFactors(2);
    end
    optimized_params_lm = x_scaled_optimized .* params.scaleFactors;
    fit_results_lm.JD = diodeModel(data_V, optimized_params_lm, config);
    fit_results_lm.resnorm = resnorm_lm;
    relative_errors_lm = abs((fit_results_lm.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    max_err_lm = max(relative_errors_lm);



    fprintf('\n第四阶段：使用trust-region-reflective算法进行拟合...\n');
    options_tr = optimoptions('lsqnonlin', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'trust-region-reflective', ...
        'FunctionTolerance', 1e-8, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'MaxFunctionEvaluations', 8000, ...
        'FiniteDifferenceType', 'central', ...
        'FiniteDifferenceStepSize', 1e-6, ...
        'MaxIterations', 4000);
    if params.lb(2) <= 0
        fprintf('警告: Rs下界为负值或零，已调整为正值\n');
        params.lb(2) = 10;
    end
    [x_scaled_optimized_tr, resnorm_tr] = runWithMultiStart(errFun, x_scaled_optimized, params.lb ./ params.scaleFactors, params.ub ./ params.scaleFactors, options_tr, optcfg);
    residual_tr = []; exitflag_tr = []; output_tr = [];

    optimized_params_tr = x_scaled_optimized_tr .* params.scaleFactors;
    fit_results_tr.JD = diodeModel(data_V, optimized_params_tr, config);
    fit_results_tr.resnorm = resnorm_tr;
    relative_errors_tr = abs((fit_results_tr.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    mean_err_tr = mean(relative_errors_tr);
    max_err_tr = max(relative_errors_tr);

    
    fprintf('\n比较两种算法的结果：\n');
    fprintf('Levenberg-Marquardt: 平均相对误差 = %.2f%%, 最大相对误差 = %.2f%%\n', ...
            mean(relative_errors_lm), max(relative_errors_lm));
    fprintf('Trust-Region-Reflective: 平均相对误差 = %.2f%%, 最大相对误差 = %.2f%%\n', ...
            mean(relative_errors_tr), max(relative_errors_tr));
            
    neg_idx = find(data_V < 0);
    pos_idx = find(data_V > 0);
    neg_err_lm = relative_errors_lm(neg_idx);
    pos_err_lm = relative_errors_lm(pos_idx);
    neg_err_tr = relative_errors_tr(neg_idx);
    pos_err_tr = relative_errors_tr(pos_idx);
    fprintf('\nLevenberg-Marquardt: 负区域误差 = %.2f%%, 正区域误差 = %.2f%%\n', mean(neg_err_lm), mean(pos_err_lm));
    fprintf('Trust-Region-Reflective: 负区域误差 = %.2f%%, 正区域误差 = %.2f%%\n', mean(neg_err_tr), mean(pos_err_tr));
    if mean(relative_errors_lm) < mean(relative_errors_tr)
        fprintf('整体表现更好: Levenberg-Marquardt算法\n');
        optimized_params = optimized_params_lm;
        fit_results.JD = fit_results_lm.JD;
        fit_results.resnorm = resnorm_lm;
        fit_results.residual = residual_lm;
        fit_results.exitflag = exitflag_lm;
        fit_results.output = output_lm;
        relative_errors = relative_errors_lm;
    else
        fprintf('整体表现更好: Trust-Region-Reflective算法\n');
        optimized_params = optimized_params_tr;
        fit_results.JD = fit_results_tr.JD;
        fit_results.resnorm = resnorm_tr;
        fit_results.residual = residual_tr;
        fit_results.exitflag = exitflag_tr;
        fit_results.output = output_tr;
        relative_errors = relative_errors_tr;
    end

    % 使用当前结果作为起点再次局部优化，以减小局部最小值的风险
    errFun_final = @(x) errorFunction(x, data_V, data_JD, params, config, config.regularization.prior);
    x_scaled_temp = optimized_params ./ params.scaleFactors;
    [x_scaled_temp, resnorm_final] = lsqnonlin(errFun_final, x_scaled_temp, params.lb ./ params.scaleFactors, params.ub ./ params.scaleFactors, options_lm);
    optimized_params = x_scaled_temp .* params.scaleFactors;
    fit_results.JD = diodeModel(data_V, optimized_params, config);
    fit_results.resnorm = resnorm_final;
    relative_errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;


    neg_errors = relative_errors(neg_idx);
    pos_errors = relative_errors(pos_idx);
    if mean(pos_errors) > 2*mean(neg_errors) && mean(pos_errors) > 10
        fprintf('\n正区域拟合效果较差，尝试单独优化正区域参数...\n');
        param_mask = [true, true, false, false, false];
        pos_errFun = @(x_opt) errorFunctionEnhancedPositive(x_opt, optimized_params ./ params.scaleFactors, param_mask, data_V, data_JD, params, config, config.regularization.prior, relative_errors);
        x0_pos_opt = optimized_params(param_mask) ./ params.scaleFactors(param_mask);
        lb_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        options_pos = optimoptions('lsqnonlin', 'Display', 'iter-detailed', 'Algorithm', 'levenberg-marquardt', 'FunctionTolerance', 1e-8, 'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, 'MaxFunctionEvaluations', 3000, 'MaxIterations', 2000);
        optcfg_enh = optcfg;
        [x_pos_opt, ~] = runWithMultiStart(pos_errFun, x0_pos_opt, lb_pos, ub_pos, options_pos, optcfg_enh);

        
        x_scaled_enhanced = optimized_params ./ params.scaleFactors;
        x_scaled_enhanced(param_mask) = x_pos_opt;
        optimized_params_enhanced = x_scaled_enhanced .* params.scaleFactors;
        fit_results_enhanced.JD = diodeModel(data_V, optimized_params_enhanced, config);
        relative_errors_enhanced = abs((fit_results_enhanced.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
        pos_errors_enhanced = relative_errors_enhanced(pos_idx);
        neg_errors_enhanced = relative_errors_enhanced(neg_idx);
        fprintf('增强前: 正区域误差 = %.2f%%, 负区域误差 = %.2f%%\n', mean(pos_errors), mean(neg_errors));
        fprintf('增强后: 正区域误差 = %.2f%%, 负区域误差 = %.2f%%\n', mean(pos_errors_enhanced), mean(neg_errors_enhanced));
        if mean(pos_errors_enhanced) < mean(pos_errors) && mean(neg_errors_enhanced) < 2*mean(neg_errors)
            fprintf('采用增强优化的结果\n');
            optimized_params = optimized_params_enhanced;
            fit_results.JD = fit_results_enhanced.JD;
            relative_errors = relative_errors_enhanced;
        else
            fprintf('保留原始优化结果\n');
        end
end

    % 计算各分量电流，用于后续分析
    currents = calculateCurrents(data_V, optimized_params, config);


    neg_idx = find(data_V < -0.1);
    neg_errors = relative_errors(neg_idx);
    fprintf('\n每个点的相对误差统计：\n');
    max_rel = max(relative_errors);
    fprintf('最大相对误差: %.2f%%\n', max_rel);
    avg_rel = mean(relative_errors);
    fprintf('平均相对误差: %.2f%%\n', avg_rel);
    prev_avg_rel = avg_rel;

    while attempt < config.optimization.max_attempts
        fprintf('误差未满足阈值，尝试再次优化 (%d/%d)...\n', ...
            attempt + 1, config.optimization.max_attempts);
        options_lm.MaxIterations = options_lm.MaxIterations * 2;
        options_lm.MaxFunctionEvaluations = options_lm.MaxFunctionEvaluations * 2;
        [optimized_params, fit_results] = final_optimization(data_V, data_JD, ...
        optimized_params ./ params.scaleFactors, params, config, options_lm, optcfg, attempt + 1, retry_count);
        relative_errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
        max_rel = max(relative_errors);
        avg_rel = mean(relative_errors);
        if prev_avg_rel - avg_rel < config.optimization.improvement_threshold
            fprintf('改进不足 %.2f%%，提前停止重试\n', config.optimization.improvement_threshold);
            break;
        end
        prev_avg_rel = avg_rel;
        attempt = attempt + 1;
    end
    
    fprintf('最终平均相对误差: %.2f%%\n', avg_rel);
    fprintf('最终最大相对误差: %.2f%%\n', max_rel);
    fprintf('中位相对误差: %.2f%%\n', median(relative_errors));
    fprintf('负电压区域平均相对误差: %.2f%%\n', mean(neg_errors));

    strong_neg_idx = find(data_V < -0.2 & data_V >= -0.5);
    if ~isempty(strong_neg_idx)
        strong_neg_errors = relative_errors(strong_neg_idx);
        fprintf('强负电压区域(-0.5到-0.2V)平均相对误差: %.2f%%\n', mean(strong_neg_errors));
    end

    pos_idx = find(data_V > 0);
    if ~isempty(pos_idx)
        pos_errors = relative_errors(pos_idx);
        fprintf('正电压区域平均相对误差: %.2f%%\n', mean(pos_errors));
        low_pos_idx = find(data_V > 0 & data_V <= 0.15);
        mid_pos_idx = find(data_V > 0.15 & data_V <= 0.25);
        high_pos_idx = find(data_V > 0.25);
        if ~isempty(low_pos_idx)
            low_pos_errors = relative_errors(low_pos_idx);
            fprintf('低正电压区域(0-0.15V)平均相对误差: %.2f%%\n', mean(low_pos_errors));
        end
        if ~isempty(mid_pos_idx)
            mid_pos_errors = relative_errors(mid_pos_idx);
            fprintf('中正电压区域(0.15-0.25V)平均相对误差: %.2f%%\n', mean(mid_pos_errors));
        end
        
        if ~isempty(high_pos_idx)
            high_pos_errors = relative_errors(high_pos_idx);
            fprintf('高正电压区域(>0.25V)平均相对误差: %.2f%%\n', mean(high_pos_errors));
        end
    end


% 针对负电压区域的误差函数
    [max_error, max_error_idx] = max(relative_errors);
    fprintf('\n误差最大的点：\n');
    fprintf('电压: %.3f V\n', data_V(max_error_idx));
    fprintf('测量电流: %.3e A\n', data_JD(max_error_idx));
    fprintf('拟合电流: %.3e A\n', fit_results.JD(max_error_idx));
    fprintf('相对误差: %.2f%%\n', max_error);

    fprintf('\n拟合参数：\n');
    fprintf('J01 = %.6e A\n', optimized_params(1));
    fprintf('Rs  = %.6e Ohm\n', optimized_params(2));
    fprintf('k   = %.6e\n', optimized_params(4));
    fprintf('J02 = %.6e A\n', optimized_params(5));
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

% 部分参数优化的误差函数
function err = errorFunctionPartial(x_opt, x0, param_mask, data_V, data_JD, params, config, prior, prev_errors)
    % 构建完整参数向量
    x_full = x0;
    x_full(param_mask) = x_opt;
    if nargin < 9
        prev_errors = [];
    end
    % 反缩放参数
    x_actual = x_full .* params.scaleFactors;
    
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
        

        % 根据上一阶段的误差加权，排除零电压附近的点
        if nargin >= 9 && ~isempty(prev_errors) && abs(data_V(i)) > 0.05
            w = 1 + abs(prev_errors(i));
            err(i) = err(i) * w;
        elseif abs(data_V(i)) > 0.05
            w = 1 + abs(err(i));
            err(i) = err(i) * w;
        end
    end

    if nargin < 7 || isempty(prior)
        prior = config.regularization.prior;
    end
    if isfield(config, 'regularization') && config.regularization.lambda > 0
        penalty = sqrt(config.regularization.lambda) * (x_actual(:) - prior(:));
        err = [err; penalty];
    end
end

% 增强正电压区域拟合的误差函数
function err = errorFunctionEnhancedPositive(x_opt, x0, param_mask, data_V, data_JD, params, config, prior, prev_errors)
    % 构建完整参数向量
    x_full = x0;
    x_full(param_mask) = x_opt;
    if nargin < 9
        prev_errors = [];
    end
    % 反缩放参数
    x_actual = x_full .* params.scaleFactors;
    
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
        
        % 根据上一阶段的误差加权，排除零电压附近的点
        if nargin >= 9 && ~isempty(prev_errors) && abs(data_V(i)) > 0.05
            w = 1 + abs(prev_errors(i));
            err(i) = err(i) * w;
        elseif abs(data_V(i)) > 0.05
            w = 1 + abs(err(i));
            err(i) = err(i) * w;
        end
    end

    if nargin < 7 || isempty(prior)
        prior = config.regularization.prior;
    end
    if isfield(config, 'regularization') && config.regularization.lambda > 0
        penalty = sqrt(config.regularization.lambda) * (x_actual(:) - prior(:));
        err = [err; penalty];
    end
end
