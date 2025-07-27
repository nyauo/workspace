function [fit_results, currents] = finalizeAdjustment(V, JD, params, config, handles)
%FINALIZEADJUSTMENT Save adjusted parameters and close figures.
%   [FIT_RESULTS, CURRENTS] = FINALIZEADJUSTMENT(...) computes the final fit,
%   saves the results and closes any open figures provided in HANDLES.

    currents = calculateCurrents(V, params, config);
    fit_results.JD = currents.total;
    fit_results.resnorm = sum(((fit_results.JD - JD) ./ (abs(JD)+eps)).^2);
    saveResults(V, JD, params, fit_results, currents);
    saveAdjustedParameters(params);

    if nargin >= 5 && isstruct(handles)
        if isfield(handles,'fitFig') && isvalid(handles.fitFig)
            close(handles.fitFig);
        end
        if isfield(handles,'errFig') && isvalid(handles.errFig)
            close(handles.errFig);
        end
    end
end
