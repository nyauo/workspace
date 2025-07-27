function updateAdjustmentPlots(handles, errors, nz_idx, currents, params, factor, avg_err)
%UPDATEADJUSTMENTPLOTS Refresh adjustment figures after a parameter change.

    set(0,'CurrentFigure',handles.fitFig);
    set(handles.h_fit,      'YData', abs(currents.total));
    set(handles.h_diode,    'YData', abs(currents.diode));
    set(handles.h_ohmic,    'YData', abs(currents.ohmic));
    set(handles.h_nonohmic, 'YData', abs(currents.nonohmic));
    xlim([-0.5 0.3]);
    ylim([1e-11 1e-3]);
    axis square;
    delete(findall(handles.fitFig,'Type','annotation'));
    annotation('textbox',[0.01,0.01,0.98,0.08], ...
        'String',sprintf('J0: %.2e  Rs: %.2e  Rsh: %.2e  k: %.2e  step: %.2f', ...
        params(1),params(2),params(3),params(4),factor), ...
        'EdgeColor','none','FontSize',10,'HorizontalAlignment','center');

    set(0,'CurrentFigure',handles.errFig);
    set(handles.h_error,'YData',errors(nz_idx));
    title(sprintf('Relative error (avg: %.2f%%)', avg_err));
    xlim([-0.5 0.3]);
    drawnow;
end
    nz_idx = V ~= 0;
    handles.h_error = bar(V(nz_idx), errors(nz_idx),'FaceColor',c_data,'EdgeColor','none');
    xlabel('Voltage (V)');
    ylabel('Relative Error (%)');
    title('');
    xlim([-0.5 0.3]);
    grid on;
end
