function [optimized_params, fit_results] = performFitting(data_V, data_JD, params, config)
    try
        % 缩放初始参数
        x0_scaled = params.x0 ./ params.scaleFactors;
        
        % 打印初始参数
        fprintf('\n初始参数（使用Lambert W函数估计）：\n');
        fprintf('J0 = %.6e A\n', params.x0(1));
        fprintf('Rs = %.6e Ohm\n', params.x0(2));
        fprintf('Rsh = %.6e Ohm\n', params.x0(3));
        fprintf('k = %.6e\n', params.x0(4));
        
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
        x0_scaled = fit_negative_region(data_V, data_JD, x0_scaled, params, config, options_lm);
        x0_scaled = fit_positive_region(data_V, data_JD, x0_scaled, params, config, options_lm);
        [optimized_params, fit_results] = final_optimization(data_V, data_JD, x0_scaled, params, config, options_lm);
    catch ME
        error('拟合过程出错: %s', ME.message);
    end
end

function x0_scaled = fit_negative_region(data_V, data_JD, x0_scaled, params, config, options_lm)
    fprintf('\n第一阶段：优化Rsh和非欧姆系数k...\n');
    neg_idx = find(data_V < -0.1);
    if ~isempty(neg_idx)
        neg_V = data_V(neg_idx);
        neg_JD = data_JD(neg_idx);
        x0_limited = x0_scaled;
        param_mask = [false, false, true, true];
        neg_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_limited, param_mask, neg_V, neg_JD, params, config);
        x0_neg_opt = x0_scaled(param_mask);
        lb_neg = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_neg = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        [x_neg_opt, ~, ~, ~, ~] = lsqnonlin(neg_errFun, x0_neg_opt, lb_neg, ub_neg, options_lm);
        x0_scaled(param_mask) = x_neg_opt;
        x_actual_neg = x0_scaled .* params.scaleFactors;
        fit_JD_neg = diodeModel(neg_V, x_actual_neg, config);
        neg_rel_errors = abs((fit_JD_neg - neg_JD) ./ (abs(neg_JD) + eps)) * 100;
        fprintf('负电压区域拟合：平均相对误差 = %.2f%%\n', mean(neg_rel_errors));
        fprintf('优化后的Rsh = %.6e Ohm\n', x_actual_neg(3));
        fprintf('优化后的k = %.6e\n', x_actual_neg(4));
    end
end

function x0_scaled = fit_positive_region(data_V, data_JD, x0_scaled, params, config, options_lm)
    fprintf('\n第二阶段：细分正电压区域优化...\n');
    low_pos_idx = find(data_V > 0 & data_V <= 0.15);
    high_pos_idx = find(data_V > 0.15);

    if ~isempty(low_pos_idx)
        fprintf('优化低正电压区域 (0 - 0.15V)...\n');
        low_pos_V = data_V(low_pos_idx);
        low_pos_JD = data_JD(low_pos_idx);
        param_mask = [true, false, false, false];
        low_pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, low_pos_V, low_pos_JD, params, config);
        x0_low_pos_opt = x0_scaled(param_mask);
        lb_low_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_low_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        [x_low_pos_opt, ~, ~, ~, ~] = lsqnonlin(low_pos_errFun, x0_low_pos_opt, lb_low_pos, ub_low_pos, options_lm);
        x0_scaled(param_mask) = x_low_pos_opt;
        x_actual_low_pos = x0_scaled .* params.scaleFactors;
        fit_JD_low_pos = diodeModel(low_pos_V, x_actual_low_pos, config);
        low_pos_rel_errors = abs((fit_JD_low_pos - low_pos_JD) ./ (abs(low_pos_JD) + eps)) * 100;
        fprintf('低正电压区域拟合：平均相对误差 = %.2f%%\n', mean(low_pos_rel_errors));
        fprintf('优化后的J0 = %.6e A\n', x_actual_low_pos(1));
    end

    if ~isempty(high_pos_idx)
        fprintf('优化高正电压区域 (>0.15V)...\n');
        high_pos_V = data_V(high_pos_idx);
        high_pos_JD = data_JD(high_pos_idx);
        param_mask = [false, true, false, false];
        high_pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, high_pos_V, high_pos_JD, params, config);
        x0_high_pos_opt = x0_scaled(param_mask);
        lb_high_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_high_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        [x_high_pos_opt, ~, ~, ~, ~] = lsqnonlin(high_pos_errFun, x0_high_pos_opt, lb_high_pos, ub_high_pos, options_lm);
        x0_scaled(param_mask) = x_high_pos_opt;
        if x0_scaled(2) * params.scaleFactors(2) <= 0
            fprintf('警告: Rs为负值或零，正在调整为正值\n');
            x0_scaled(2) = params.lb(2) / params.scaleFactors(2);
        end
        
        x_actual_high_pos = x0_scaled .* params.scaleFactors;
        fit_JD_high_pos = diodeModel(high_pos_V, x_actual_high_pos, config);
        high_pos_rel_errors = abs((fit_JD_high_pos - high_pos_JD) ./ (abs(high_pos_JD) + eps)) * 100;
        fprintf('高正电压区域拟合：平均相对误差 = %.2f%%\n', mean(high_pos_rel_errors));
        fprintf('优化后的Rs = %.6e Ohm\n', x_actual_high_pos(2));
    end

    pos_idx = find(data_V > 0);
    if ~isempty(pos_idx)
        fprintf('\n综合优化正电压区域...\n');
        pos_V = data_V(pos_idx);
        pos_JD = data_JD(pos_idx);
        param_mask = [true, true, false, false];
        pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, pos_V, pos_JD, params, config);
        x0_pos_opt = x0_scaled(param_mask);
        lb_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        [x_pos_opt, ~, ~, ~, ~] = lsqnonlin(pos_errFun, x0_pos_opt, lb_pos, ub_pos, options_lm);
        x0_scaled(param_mask) = x_pos_opt;
        if x0_scaled(2) * params.scaleFactors(2) <= 0
            fprintf('警告: Rs为负值或零，正在调整为正值\n');
            x0_scaled(2) = params.lb(2) / params.scaleFactors(2);
        end
        
        x_actual_pos = x0_scaled .* params.scaleFactors;
        fit_JD_pos = diodeModel(pos_V, x_actual_pos, config);
        pos_rel_errors = abs((fit_JD_pos - pos_JD) ./ (abs(pos_JD) + eps)) * 100;
        fprintf('正电压区域综合拟合：平均相对误差 = %.2f%%\n', mean(pos_rel_errors));
        fprintf('优化后的J0 = %.6e A\n', x_actual_pos(1));
        fprintf('优化后的Rs = %.6e Ohm\n', x_actual_pos(2));
    end
end

function [optimized_params, fit_results] = final_optimization(data_V, data_JD, x0_scaled, params, config, options_lm)
    fprintf('\n第三阶段：全区域拟合...\n');
    errFun = @(x) errorFunction(x, data_V, data_JD, params, config);
    [x_scaled_optimized, resnorm_lm, residual_lm, exitflag_lm, output_lm] = lsqnonlin(errFun, x0_scaled, [], [], options_lm);
    if x_scaled_optimized(2) * params.scaleFactors(2) <= 0
        fprintf('警告: LM算法产生了负值或零的Rs，正在调整为正值\n');
        x_scaled_optimized(2) = params.lb(2) / params.scaleFactors(2);
    end
    optimized_params_lm = x_scaled_optimized .* params.scaleFactors;
    fit_results_lm.JD = diodeModel(data_V, optimized_params_lm, config);
    fit_results_lm.resnorm = resnorm_lm;
    relative_errors_lm = abs((fit_results_lm.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;

    fprintf('\n第四阶段：使用trust-region-reflective算法进行拟合...\n');
    options_tr = optimoptions('lsqnonlin', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'trust-region-reflective', ...
        'FunctionTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'MaxFunctionEvaluations', 8000, ...
        'FiniteDifferenceType', 'central', ...
        'FiniteDifferenceStepSize', 1e-6, ...
        'MaxIterations', 4000);
    if params.lb(2) <= 0
        fprintf('警告: Rs下界为负值或零，已调整为正值\n');
        params.lb(2) = 10;
    end
    [x_scaled_optimized_tr, resnorm_tr, residual_tr, exitflag_tr, output_tr] = ...
        lsqnonlin(errFun, x_scaled_optimized, params.lb ./ params.scaleFactors, params.ub ./ params.scaleFactors, options_tr);
    if x_scaled_optimized_tr(2) * params.scaleFactors(2) <= 0
        fprintf('警告: TR算法产生了负值或零的Rs，正在调整为正值\n');
        x_scaled_optimized_tr(2) = params.lb(2) / params.scaleFactors(2);
    end
    optimized_params_tr = x_scaled_optimized_tr .* params.scaleFactors;
    fit_results_tr.JD = diodeModel(data_V, optimized_params_tr, config);
    fit_results_tr.resnorm = resnorm_tr;
    relative_errors_tr = abs((fit_results_tr.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;

    fprintf('\n比较两种算法的结果：\n');
    fprintf('Levenberg-Marquardt: 平均相对误差 = %.2f%%\n', mean(relative_errors_lm));
    fprintf('Trust-Region-Reflective: 平均相对误差 = %.2f%%\n', mean(relative_errors_tr));

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

    if optimized_params(2) <= 0
        fprintf('警告: 最终拟合结果中Rs为负值或零，已调整为正值\n');
        optimized_params(2) = max(params.lb(2), 10);
        fit_results.JD = diodeModel(data_V, optimized_params, config);
        relative_errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    end

    neg_errors = relative_errors(neg_idx);
    pos_errors = relative_errors(pos_idx);
    if mean(pos_errors) > 2*mean(neg_errors) && mean(pos_errors) > 10
        fprintf('\n正区域拟合效果较差，尝试单独优化正区域参数...\n');
        param_mask = [true, true, false, false];
        pos_errFun = @(x_opt) errorFunctionEnhancedPositive(x_opt, optimized_params ./ params.scaleFactors, param_mask, data_V, data_JD, params, config);
        x0_pos_opt = optimized_params(param_mask) ./ params.scaleFactors(param_mask);
        lb_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        options_pos = optimoptions('lsqnonlin', 'Display', 'iter-detailed', 'Algorithm', 'levenberg-marquardt', 'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'MaxFunctionEvaluations', 3000, 'MaxIterations', 2000);
        [x_pos_opt, ~, ~, ~, ~] = lsqnonlin(pos_errFun, x0_pos_opt, [], [], options_pos);
        if x_pos_opt(2) * params.scaleFactors(2) <= 0
            fprintf('警告: 正区域优化产生了负值或零的Rs，正在调整为正值\n');
            x_pos_opt(2) = lb_pos(2);
        end
        
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

    % 计算欧姆与非欧姆电流的平均绝对值并检查其比值
    % 如果非欧姆电流显著大于欧姆电流，可能表明模型不符合物理约束
    currents = calculateCurrents(data_V, optimized_params, config);
    mean_ohmic = mean(abs(currents.ohmic));
    mean_nonohmic = mean(abs(currents.nonohmic));
    ratio = mean_nonohmic / (mean_ohmic + eps);
    ratio_threshold = 1e3;
    if ratio > ratio_threshold
        fprintf(['警告: 非欧姆电流与欧姆电流的平均值比率为 %.2e, 超过阈值 %.0f, ' ...
                '尝试在收紧k上界后重新拟合...\n'], ratio, ratio_threshold);
        params.ub(4) = min(params.ub(4), optimized_params(4));
        [optimized_params, fit_results] = final_optimization(data_V, data_JD, ...
            optimized_params ./ params.scaleFactors, params, config, options_lm);
        relative_errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    end

    neg_idx = find(data_V < -0.1);
    neg_errors = relative_errors(neg_idx);
    fprintf('\n每个点的相对误差统计：\n');
    fprintf('最大相对误差: %.2f%%\n', max(relative_errors));
    fprintf('平均相对误差: %.2f%%\n', mean(relative_errors));
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
    fprintf('J0 = %.6e A\n', optimized_params(1));
    fprintf('Rs = %.6e Ohm\n', optimized_params(2));
    fprintf('Rsh = %.6e Ohm\n', optimized_params(3));
    fprintf('k = %.6e\n', optimized_params(4));
end
function err = errorFunctionNegative(x, data_V, data_JD, params, config)
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
end

% 部分参数优化的误差函数
function err = errorFunctionPartial(x_opt, x0, param_mask, data_V, data_JD, params, config)
    % 构建完整参数向量
    x_full = x0;
    x_full(param_mask) = x_opt;
    
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
        
        % 针对不同电压区域应用不同权重
        if data_V(i) < -0.3
            err(i) = err(i) * 5; % 强负电压区域
        elseif data_V(i) < -0.1
            err(i) = err(i) * 3; % 弱负电压区域
        elseif data_V(i) < 0.1
            err(i) = err(i) * 2; % 零点附近
        end
    end
end

% 增强正电压区域拟合的误差函数
function err = errorFunctionEnhancedPositive(x_opt, x0, param_mask, data_V, data_JD, params, config)
    % 构建完整参数向量
    x_full = x0;
    x_full(param_mask) = x_opt;
    
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
        
        % 特别加强正电压区域权重，特别是拟合不佳的高电压区域
        if data_V(i) > 0.25
            err(i) = err(i) * 8;  % 特别强调高正电压区域
        elseif data_V(i) > 0.15
            err(i) = err(i) * 5;  % 中正电压区域
        elseif data_V(i) > 0
            err(i) = err(i) * 3;  % 低正电压区域
        elseif data_V(i) < -0.3
            err(i) = err(i) * 0.5; % 弱化强负电压区域的影响
        end
    end
end
