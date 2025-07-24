function err = errorFunction(x, data_V, data_JD, params, config, prior)
    % 反缩放参数
    x_actual = x .* params.scaleFactors;
    
    % 计算模型预测值
    predicted = diodeModel(data_V, x_actual, config);
    
    % 计算误差，使用对数空间误差
    % 这对于跨多个数量级的数据特别有效
    err = zeros(size(data_JD));
    
    for i = 1:length(data_JD)
        % 确保数据为正值以便取对数（对于二极管电流，大多数时候是正的）
        actual_abs = abs(data_JD(i));
        pred_abs = abs(predicted(i));
        
        % 阈值，避免对近零值取对数
        threshold = 1e-12;
        
        % 处理小于阈值的情况
        if actual_abs < threshold || pred_abs < threshold
            % 对小值直接使用标准化误差
            err(i) = (predicted(i) - data_JD(i)) / max(1e-12, max(max(abs(data_JD))));
        else
            % 计算对数空间的误差
            log_actual = log10(actual_abs);
            log_pred = log10(pred_abs);
            err(i) = log_pred - log_actual;
            
            % 保持符号一致性
            if sign(predicted(i)) ~= sign(data_JD(i)) && abs(predicted(i)) > threshold && abs(data_JD(i)) > threshold
                err(i) = err(i) * 3; % 对符号不同的情况增加误差权重
            end
        end
          % 根据误差大小动态加权，零电压附近不加权
        if abs(data_V(i)) > 0.05
            w = 1 + abs(err(i));
            err(i) = err(i) * w;
        end
    end      

    % Append L2 penalty based on parameter priors if enabled
    if nargin < 6 || isempty(prior)
        prior = config.regularization.prior;
    end
    if isfield(config, 'regularization') && config.regularization.lambda > 0
        penalty = sqrt(config.regularization.lambda) * (x_actual(:) - prior(:));
        err = [err; penalty];
    end
end
