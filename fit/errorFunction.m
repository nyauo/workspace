function err = errorFunction(x, data_V, data_JD, params, config, prior, param_mask, x0, prev_errors, sign_weight)
%ERRORFUNCTION Unified error computation for diode model fitting.
%   ERR = ERRORFUNCTION(X, V, JD, PARAMS, CONFIG) computes the residual
%   vector between the measured current JD and the model prediction using
%   parameters X.  Optional arguments allow partial optimisation and
%   custom weighting:
%       PARAM_MASK - logical vector indicating which parameters of X0 are
%                    replaced by X when performing a partial update.
%       X0         - full parameter vector used together with PARAM_MASK.
%       PREV_ERRORS- previous residuals used for adaptive weighting.
%       SIGN_WEIGHT- factor applied when the predicted and measured
%                    currents have different signs (default 3 for full
%                    optimisation, 4 when PARAM_MASK is supplied).

    if nargin < 6 || isempty(prior)
        prior = config.regularization.prior;
    end
    if nargin < 7, param_mask = []; end
    if nargin < 8, x0 = []; end
    if nargin < 9, prev_errors = []; end
    if nargin < 10 || isempty(sign_weight)
        if ~isempty(param_mask)
            sign_weight = 4; % behaviour of previous partial variants
        else
            sign_weight = 3; % default for full optimisation
        end
    end

    if ~isempty(param_mask)
        x_full = x0;
        x_full(param_mask) = x;
    else
        x_full = x;
    end

    % 反缩放参数
    x_actual = x_full .* params.scaleFactors;

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

    if ~isempty(param_mask)
        sign_mismatch = sign(predicted) ~= sign(data_JD);
        err(sign_mismatch & mask_log) = err(sign_mismatch & mask_log) * sign_weight;
        idx = abs(data_V) > 0.05;
        if ~isempty(prev_errors)
            err(idx) = err(idx) .* (1 + abs(prev_errors(idx)));
        else
            err(idx) = err(idx) .* (1 + abs(err(idx)));
        end
    else
        sign_mismatch = sign(predicted) ~= sign(data_JD) & pred_abs > threshold & actual_abs > threshold;
        err(sign_mismatch) = err(sign_mismatch) * sign_weight;
        idx = abs(data_V) > 0.05;
        err(idx) = err(idx) .* (1 + abs(err(idx)));
    end

    if isfield(config, 'regularization') && config.regularization.lambda > 0
        penalty = sqrt(config.regularization.lambda) * (x_actual(:) - prior(:));
        err = [err; penalty];
    end
end
