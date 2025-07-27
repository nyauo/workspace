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
    threshold = 1e-12;
    max_val = max(1e-12, max(abs(data_JD(:))));
    actual_abs = abs(data_JD);
    pred_abs = abs(predicted);
    err = zeros(size(data_JD));
    mask_small = actual_abs < threshold | pred_abs < threshold;
    err(mask_small) = (predicted(mask_small) - data_JD(mask_small)) / max_val;

    mask_log = ~mask_small;
    err(mask_log) = log10(pred_abs(mask_log)) - log10(actual_abs(mask_log));

    sign_mismatch = sign(predicted) ~= sign(data_JD);
    err(sign_mismatch & mask_log) = err(sign_mismatch & mask_log) * 4;

    idx = abs(data_V) > 0.05;
    if nargin >= 9 && ~isempty(prev_errors)
        err(idx) = err(idx) .* (1 + abs(prev_errors(idx)));
    else
        err(idx) = err(idx) .* (1 + abs(err(idx)));
    end

    if nargin < 7 || isempty(prior)
        prior = config.regularization.prior;
    end
    if isfield(config, 'regularization') && config.regularization.lambda > 0
        penalty = sqrt(config.regularization.lambda) * (x_actual(:) - prior(:));
        err = [err; penalty];
    end
end
