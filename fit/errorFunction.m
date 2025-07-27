function err = errorFunction(x, data_V, data_JD, params, config, prior)
    % 反缩放参数
    x_actual = x .* params.scaleFactors;
    
    % 计算模型预测值
    predicted = diodeModel(data_V, x_actual, config);
        
    threshold = 1e-12;
    max_val = max(1e-12, max(abs(data_JD(:))));
    actual_abs = abs(data_JD);
    pred_abs = abs(predicted);

    % 对非常小的值直接使用标准化误差
    mask_small = actual_abs < threshold | pred_abs < threshold;
    err = zeros(size(data_JD));
    err(mask_small) = (predicted(mask_small) - data_JD(mask_small)) / max_val;

    % 计算对数空间的误差
    mask_log = ~mask_small;
    err(mask_log) = log10(pred_abs(mask_log)) - log10(actual_abs(mask_log));

    % 保持符号一致性，对符号不同的情况增加误差权重
    sign_mismatch = sign(predicted) ~= sign(data_JD) & pred_abs > threshold & actual_abs > threshold;
    err(sign_mismatch) = err(sign_mismatch) * 3;

    % 根据误差大小动态加权，零电压附近不加权
    idx = abs(data_V) > 0.05;
    err(idx) = err(idx) .* (1 + abs(err(idx)));
    % 计算误差，使用对数空间误差
    % 这对于跨多个数量级的数据特别有效
   

    % Append L2 penalty based on parameter priors if enabled
    if nargin < 6 || isempty(prior)
        prior = config.regularization.prior;
    end
    if isfield(config, 'regularization') && config.regularization.lambda > 0
        penalty = sqrt(config.regularization.lambda) * (x_actual(:) - prior(:));
        err = [err; penalty];
    end
end
