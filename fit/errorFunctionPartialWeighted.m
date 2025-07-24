function err = errorFunctionPartialWeighted(x_opt, x0, param_mask, data_V, data_JD, params, config, prior, prev_errors)
    % 构建完整参数向量
    if nargin < 9
        prev_errors = [];
    end
    x_full = x0;
    x_full(param_mask) = x_opt;
    x_actual = x_full .* params.scaleFactors;

    % 计算模型预测值
    predicted = diodeModel(data_V, x_actual, config);

    err = zeros(size(data_JD));
    for i = 1:length(data_JD)
        actual_abs = abs(data_JD(i));
        pred_abs = abs(predicted(i));
        threshold = 1e-12;
        if actual_abs < threshold || pred_abs < threshold
            err(i) = (predicted(i) - data_JD(i)) / max(1e-12, max(max(abs(data_JD))));
        else
            log_actual = log10(actual_abs);
            log_pred = log10(pred_abs);
            err(i) = log_pred - log_actual;
            if sign(predicted(i)) ~= sign(data_JD(i))
                err(i) = err(i) * 4;
            end
        end
        if ~isempty(prev_errors) && abs(data_V(i)) > 0.05
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
