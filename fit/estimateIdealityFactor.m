function n_est = estimateIdealityFactor(V, JD, config)
%ESTIMATEIDEALITYFACTOR Estimate diode ideality factor from measured data
%   n_est = ESTIMATEIDEALITYFACTOR(V, JD, config) computes the local slope of
%   log(abs(JD)) with respect to V for V > 0 and converts it to an estimate of
%   the diode ideality factor. The returned value is the median of the
%   point-wise estimates clipped to the range [1, 2].

    % Select positive voltage region
    idx = find(V > 0);
    if numel(idx) < 2
        % Not enough points - fall back to configured value
        n_est = config.physics.n;
        return;
    end

    V_pos = V(idx);
    JD_pos = JD(idx);

    % Derivative of log(abs(JD)) w.r.t voltage
    slopes = diff(log(abs(JD_pos))) ./ diff(V_pos);

    % Convert slopes to ideality factor estimates
    n_vals = config.physics.q ./ (config.physics.kb * config.physics.T .* slopes);

    % Median estimate clipped to [1, 2]
    n_est = median(n_vals);
    n_est = max(1, min(2, n_est));
end
