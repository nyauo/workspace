function JD = diodeModel(V, params, config)    
    % 检查Rs参数是否为正值
    if params(2) <= 0
        error('物理参数错误: Rs必须为正值 (当前值: %.6e)', params(2));
    end
    
    % 设置fsolve选项
    options = optimoptions('fsolve', ...
        'Display', 'off', ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'StepTolerance', 1e-12);
    

    % 计算初始猜测（向量形式）
    initial_guess = zeros(size(V));
    neg_idx = V < 0;
    pos_idx = ~neg_idx;

    initial_guess(neg_idx) = (V(neg_idx) ./ params(3)) + ...
        params(4) * (abs(V(neg_idx)).^config.physics.m) .* sign(V(neg_idx));
    initial_guess(pos_idx) = params(1) * (exp(config.physics.A * V(pos_idx) / config.physics.n) - 1);

    % 向量化残差函数
    residual = @(J) params(1) * (exp(config.physics.A * (V - J .* params(2)) / config.physics.n) - 1) + ...
                    (V - J .* params(2)) ./ params(3) + ...
                    params(4) * (abs(V - J .* params(2)).^config.physics.m) .* sign(V - J .* params(2)) - J;

    % 求解所有电压点的电流
    JD = fsolve(residual, initial_guess, options);
end
