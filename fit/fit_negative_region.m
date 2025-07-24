function [x0_scaled, rel_errors] = fit_negative_region(data_V, data_JD, x0_scaled, params, config, options_lm, optcfg)
    fprintf('\n第一阶段：优化Rsh和非欧姆系数k...\n');
    neg_idx = find(data_V < -0.1);
    if ~isempty(neg_idx)
        neg_V = data_V(neg_idx);
        neg_JD = data_JD(neg_idx);
        x0_limited = x0_scaled;
        param_mask = [false, false, true, true, false];
        neg_errFun = @(x_opt) errorFunctionPartialWeighted(x_opt, x0_limited, param_mask, neg_V, neg_JD, params, config, config.regularization.prior, []);
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
