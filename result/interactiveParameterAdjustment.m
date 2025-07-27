function [adjusted_params, fit_results] = interactiveParameterAdjustment(V, JD, initial_params, config)
%INTERACTIVEPARAMETERADJUSTMENT Command-line wrapper for parameter tweaking.
%   This routine loops over user choices and calls helper functions to
%   update the adjustment figures.
% Copy initial parameters and set default adjustment step    adjusted_params = initial_params;
adjusted_params = initial_params;
adjustment_factor = 1.0;
    
    if adjusted_params(2) <= 0
        fprintf('Warning: initial Rs non-positive, setting to 10.\n');
        adjusted_params(2) = 10;
    end

    currents = calculateCurrents(V, adjusted_params, config);
    errors   = abs((currents.total - JD) ./ (abs(JD) + eps)) * 100;
    [handles, nz_idx] = setupAdjustmentFigures(V, JD, currents, errors);
    avg_error = mean(errors(nz_idx));
    updateAdjustmentPlots(handles, errors, nz_idx, currents, adjusted_params, adjustment_factor, avg_error);

    
    
    while true
        fprintf('\nCurrent: J0=%.2e  Rs=%.2e  Rsh=%.2e  k=%.2e\n', ...
            adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4));
        fprintf('Average relative error: %.2f%%\n', avg_error);
        fprintf('\n1:+J0 2:-J0 3:+Rs 4:-Rs 5:+Rsh 6:-Rsh 7:+k 8:-k\n');
        fprintf('9:change step (%.2f) 0:finish\n', adjustment_factor);
        choice = str2double(input('Choose (0-9): ','s'));
        if isnan(choice)
            fprintf('Please enter a number between 0 and 9.\n');
            continue;
        end
        
        if choice == 0
            [fit_results, ~] = finalizeAdjustment(V, JD, adjusted_params, config, handles);
            break;
        elseif choice == 9
            nf = str2double(input(sprintf('New step (current %.2f): ', adjustment_factor),'s'));
            if ~isnan(nf) && nf > 0
                adjustment_factor = nf;
            else
                fprintf('Invalid step, keeping %.2f\n', adjustment_factor);
            end
            continue;
        elseif choice >=1 && choice <=8
            [adjusted_params, currents, errors, avg_error] = applyAdjustment(choice, adjusted_params, adjustment_factor, V, JD, config, nz_idx);
            updateAdjustmentPlots(handles, errors, nz_idx, currents, adjusted_params, adjustment_factor, avg_error);
        else
            fprintf('Invalid choice.\n');
        end
    end
    
    if adjusted_params(2) <= 0
        adjusted_params(2) = 1;
    end

    if ~exist('fit_results','var')
        [fit_results, ~] = finalizeAdjustment(V, JD, adjusted_params, config, handles);
    end
end
