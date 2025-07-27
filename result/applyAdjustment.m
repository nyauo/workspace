function [params, currents, errors, avg_err] = applyAdjustment(choice, params, factor, V, JD, config, nz_idx)
%APPLYADJUSTMENT Modify parameters according to a selected action.
%   CHOICE specifies which parameter to change and in what direction.

    param_idx = ceil(choice/2);
    if mod(choice,2)==1
        direction = 1;
    else
        direction = -1;
    end

    delta = params(param_idx) * 0.1 * factor * direction;
    params(param_idx) = params(param_idx) + delta;

    if param_idx == 1
        params(param_idx) = max(1e-12, params(param_idx));
    elseif param_idx == 2
        params(param_idx) = max(1, params(param_idx));
        if params(param_idx) <= 0
            params(param_idx) = 1;
        end
    elseif param_idx == 3
        params(param_idx) = max(1e4, params(param_idx));
    elseif param_idx == 4
        params(param_idx) = max(1e-10, params(param_idx));
    end

    currents = calculateCurrents(V, params, config);
    fit_JD = currents.total;
    errors = abs((fit_JD - JD) ./ (abs(JD)+eps)) * 100;
    avg_err = mean(errors(nz_idx));
end
