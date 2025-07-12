function config = loadConfig()
    config.physics = struct(...
        'q', 1.602e-19, ...    % 电子电荷 (C)
        'kb', 1.38e-23, ...    % 波尔兹曼常数 (J/K)
        'T', 300, ...          % 温度 (K)
        'n', 1.4, ...          % Ideality factor - 增大以更好适应负区域
        'm', 2.4 ...           % Exponent for recombination - 对非欧姆项的指数调整
    );
    
    config.physics.A = config.physics.q / (config.physics.kb * config.physics.T);
    
    % 设置不同电压区域的参数
    config.fitting = struct(...
        'neg_voltage_threshold', -0.2, ... % 负电压区域阈值
        'pos_voltage_threshold', 0.1 ...   % 正电压区域阈值
    );

    % 优化相关设置
    config.optimization = struct(...
        'multistart_points', 5, ...
        'method', 'multistart', ...
        'target_rel_error', 2, ... % 收敛判据: 平均相对误差 < 2%
        'target_max_error', 5, ...  % 收敛判据: 最大相对误差 < 5%
        'max_attempts', 3 ...      % 达到阈值前最多重新优化次数
    );
    % Regularization settings for parameter priors
    config.regularization = struct(...
        'lambda', 0, ...                % penalty strength
        'prior', [0, 0, 0, 0] ...       % prior parameter values
    );    
end
