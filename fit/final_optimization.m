function [optimized_params, fit_results] = final_optimization(data_V, data_JD, x0_scaled, params, config, options_lm, optcfg, attempt)
    if nargin < 8
        attempt = 1;
    end
    
    prev_avg_rel = Inf;
    while true
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
        pos_errFun = @(x_opt) errorFunctionPartialWeighted(x_opt, optimized_params ./ params.scaleFactors, param_mask, data_V, data_JD, params, config, config.regularization.prior, relative_errors);
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
        


    neg_idx = find(data_V < -0.1);
    neg_errors = relative_errors(neg_idx);
    fprintf('\n每个点的相对误差统计：\n');
    max_rel = max(relative_errors);
    fprintf('最大相对误差: %.2f%%\n', max_rel);
    avg_rel = mean(relative_errors);
    fprintf('平均相对误差: %.2f%%\n', avg_rel);
    if avg_rel < config.optimization.termination_avg_error
        fprintf('平均误差 %.2f%% 低于阈值，结束优化。\n', avg_rel);
        break;
    end
    prev_avg_rel = avg_rel;

    while attempt < config.optimization.max_attempts
        fprintf('误差未满足阈值，尝试再次优化 (%d/%d)...\n', ...
            attempt + 1, config.optimization.max_attempts);
        options_lm.MaxIterations = options_lm.MaxIterations * 2;
        options_lm.MaxFunctionEvaluations = options_lm.MaxFunctionEvaluations * 2;
        [optimized_params, fit_results] = final_optimization(data_V, data_JD, ...
        optimized_params ./ params.scaleFactors, params, config, options_lm, optcfg, attempt + 1);
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
    
    if avg_rel < config.optimization.termination_avg_error
        fprintf('平均误差 %.2f%% 低于阈值，结束优化。\n', avg_rel);
        break;
    end

    if attempt >= config.optimization.max_attempts && avg_rel >= config.optimization.termination_avg_error
        fprintf('已达到最大尝试次数 %d，仍未达到误差阈值，停止优化。\n', config.optimization.max_attempts);
        break;
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
end
