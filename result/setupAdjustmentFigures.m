function [handles, nz_idx] = setupAdjustmentFigures(V, JD, currents, errors)
%SETUPADJUSTMENTFIGURES Create plots for parameter tweaking.
%   [HANDLES, NZ_IDX] = SETUPADJUSTMENTFIGURES(V, JD, CURRENTS, ERRORS) creates
%   the fitting and error figures used during interactive adjustment and
%   returns the handles for further updates. NZ_IDX contains the indices of
%   non-zero voltages used when computing error statistics.

    handles.fitFig = figure('Name','Adjustment Results','Position',[100 100 600 600]);
    handles.errFig = figure('Name','Error Analysis','Position',[750 100 600 400]);

    figure(handles.fitFig);
    c_data   = [107,174,214]/255; % #6BAED6
    c_total  = [251,106, 74]/255; % #FB6A4A
    c_ohmic  = [144,186, 72]/255; % #90BA48
    c_diode  = [ 19,106,238]/255; % #136AEE
    c_nonohm = [223, 66,227]/255; % #DF42E3

    handles.h_data = semilogy(V, abs(JD), 'o', 'Color', c_data, 'DisplayName','Measured');
    hold on;
    handles.h_fit  = semilogy(V, abs(currents.total), 'o', 'Color', c_total, 'DisplayName','Fitted');
    handles.h_diode = semilogy(V, abs(currents.diode),'--','Color',c_diode,'DisplayName','Diode');
    handles.h_ohmic = semilogy(V, abs(currents.ohmic),'--','Color',c_ohmic,'DisplayName','Ohmic');
    handles.h_nonohmic = semilogy(V, abs(currents.nonohmic),'--','Color',c_nonohm,'DisplayName','Non-ohmic');
    xlim([-0.5 0.3]);
    ylim([1e-11 1e-3]);
    axis square;
    xlabel('Voltage (V)');
    ylabel('Current Density (A)');
    title('I-V Characteristics (log scale)');
    legend('Location','best');
    grid on;

    figure(handles.errFig);
    nz_idx = V ~= 0;
    handles.h_error = bar(V(nz_idx), errors(nz_idx),'FaceColor',c_data,'EdgeColor','none');
    xlabel('Voltage (V)');
    ylabel('Relative Error (%)');
    title('');
    xlim([-0.5 0.3]);
    grid on;
end
