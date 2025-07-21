function config = loadConfig()
    config.physics = struct(...
        'q', 1.602e-19, ...    % 电子电荷 (C)
        'kb', 1.38e-23, ...    % 波尔兹曼常数 (J/K)
        'T', 300, ...          % 温度 (K)
        'n', 1, ...          % Ideality factor - 增大以更好适应负区域
        'n2', 2, ...           % tunnelling ideality factor
        'm', 2.4 ...           % Exponent for recombination - 对非欧姆项的指数调整
    );
    
    config.physics.A = config.physics.q / (config.physics.kb * config.physics.T);
    config.physics.A2 = config.physics.q / (config.physics.kb * config.physics.T * config.physics.n2);
    % 设置不同电压区域的参数
    config.fitting = struct(...
        'neg_voltage_threshold', -0.2, ... % 负电压区域阈值
        'pos_voltage_threshold', 0.1 ...   % 正电压区域阈值
    );

    % 优化相关设置
    config.optimization = struct(...
        'multistart_points', 5, ...
        'method', 'multistart', ...
        'max_attempts', 3, ...     % 最大优化尝试次数
        'max_retries', 2, ...      % 非欧姆比率过高时的最大递归次数
        'improvement_threshold', 0.1, ... % 最小改进幅度 (百分比)
        'm_range', [2, 3], ...     % 调整m的范围
        'm_step', 0.1 ...          % 调整m的步长
    );
    % Regularization settings for parameter priors
    config.regularization = struct(...
        'lambda', 0, ...                % penalty strength
        'prior', [0, 0, 0, 0, 0] ...    % prior parameter values
    );
    
    % Parallel computing options
    config.parallel = struct(...
        'use', false, ...     % whether to run optimisations in parallel
        'poolSize', [] ...    % optional pool size for parpool
    );
end
