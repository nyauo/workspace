function [x0_scaled, rel_errors] = fit_positive_region(data_V, data_JD, x0_scaled, params, config, options_lm, prev_errors)
    fprintf('\n第二阶段：细分正电压区域优化...\n');
    low_pos_idx = find(data_V > 0 & data_V <= 0.15);
    high_pos_idx = find(data_V > 0.15);

    if ~isempty(low_pos_idx)
        fprintf('优化低正电压区域 (0 - 0.15V)...\n');
        low_pos_V = data_V(low_pos_idx);
        low_pos_JD = data_JD(low_pos_idx);
        param_mask = [true, false, false, false];
        low_pos_prev = prev_errors(low_pos_idx);
        low_pos_errFun = @(x_opt) errorFunction(x_opt, low_pos_V, low_pos_JD, params, config, config.regularization.prior, param_mask, x0_scaled, low_pos_prev, 4);
        x0_low_pos_opt = x0_scaled(param_mask);
        lb_low_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_low_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        [x_low_pos_opt, ~] = runWithMultiStart(low_pos_errFun, x0_low_pos_opt, lb_low_pos, ub_low_pos, options_lm, config.optimization);
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
        param_mask = [false, true, false, false];
        high_pos_prev = prev_errors(high_pos_idx);
        high_pos_errFun = @(x_opt) errorFunctionPartial(x_opt, x0_scaled, param_mask, high_pos_V, high_pos_JD, params, config, config.regularization.prior, high_pos_prev, 4);
        x0_high_pos_opt = x0_scaled(param_mask);
        lb_high_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_high_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        [x_high_pos_opt, ~] = runWithMultiStart(high_pos_errFun, x0_high_pos_opt, lb_high_pos, ub_high_pos, options_lm, config.optimization);
        x0_scaled(param_mask) = x_high_pos_opt;
        % 更新误差供下一小阶段使用
        x_all = x0_scaled .* params.scaleFactors;
        fit_all = diodeModel(data_V, x_all, config);
        prev_errors = abs((fit_all - data_JD) ./ (abs(data_JD) + eps));
        if x0_scaled(2) * params.scaleFactors(2) <= 0
            fprintf('警告: Rs为负值或零，正在调整为正值\n');
            x0_scaled(2) = params.lb(2) / params.scaleFactors(2);
        end
        
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
        param_mask = [true, true, false, false];
        pos_prev = prev_errors(pos_idx);
        pos_errFun = @(x_opt) errorFunction(x_opt, pos_V, pos_JD, params, config, config.regularization.prior, param_mask, x0_scaled, pos_prev, 4);
        x0_pos_opt = x0_scaled(param_mask);
        lb_pos = params.lb(param_mask) ./ params.scaleFactors(param_mask);
        ub_pos = params.ub(param_mask) ./ params.scaleFactors(param_mask);
        [x_pos_opt, ~] = runWithMultiStart(pos_errFun, x0_pos_opt, lb_pos, ub_pos, options_lm, config.optimization);
        x0_scaled(param_mask) = x_pos_opt;
        if x0_scaled(2) * params.scaleFactors(2) <= 0
            fprintf('警告: Rs为负值或零，正在调整为正值\n');
            x0_scaled(2) = params.lb(2) / params.scaleFactors(2);
        end
        
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
