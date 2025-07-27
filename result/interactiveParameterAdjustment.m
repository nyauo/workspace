function [adjusted_params, fit_results] = interactiveParameterAdjustment(data_V, data_JD, initial_params, config)
    %INTERACTIVEPARAMETERADJUSTMENT Manually tweak model parameters.
    %   [ADJUSTED_PARAMS, FIT_RESULTS] = INTERACTIVEPARAMETERADJUSTMENT(V, JD,
    %   INITIAL_PARAMS, CONFIG) opens an interactive loop that allows the user
    %   to adjust parameters and immediately see the effect on the fit.

    % Copy initial parameters and set default adjustment step
    adjusted_params = initial_params;
    adjustment_factor = 1.0;
    
    % 确保初始参数物理合理性
    if adjusted_params(2) <= 0  % Rs必须为正
        fprintf('警告: 初始Rs为负值或零，已自动调整为正值\n');
        adjusted_params(2) = 10; % 使用一个合理的默认值
    end
    
    % 计算初始拟合和误差
    %fit_results.JD = diodeModel(data_V, adjusted_params, config);
    currents = calculateCurrents(data_V, adjusted_params, config);
    fit_results.JD = currents.total;
    errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    nz_idx = data_V ~= 0;
    avg_error = mean(errors(nz_idx));

    % Create figures
    fitFig = figure('Name','拟合结果','Position',[100 100 600 600]);
    errFig = figure('Name','误差分析','Position',[750 100 600 400]);

    % Colours
    figure(fitFig);
    c_data   = [107,174,214]/255; % #6BAED6
    c_total  = [251,106, 74]/255; % #FB6A4A
    c_ohmic  = [144,186, 72]/255; % #90BA48
    c_diode  = [ 19,106,238]/255; % #136AEE
    c_nonohm = [223, 66,227]/255; % #DF42E3

    % Semilog I-V plot
    h_data = semilogy(data_V, abs(data_JD), 'o', 'Color', c_data, 'DisplayName', '测量数据');
    hold on;
    h_fit = semilogy(data_V, abs(currents.total), 'o', 'Color', c_total, 'DisplayName', '拟合结果');
    h_diode = semilogy(data_V, abs(currents.diode), '--', 'Color', c_diode, 'DisplayName', '二极管电流');
    h_ohmic = semilogy(data_V, abs(currents.ohmic), '--', 'Color', c_ohmic, 'DisplayName', '欧姆电流');
    h_nonohmic = semilogy(data_V, abs(currents.nonohmic), '--', 'Color', c_nonohm, 'DisplayName', '非欧姆电流');
    xlim([-0.5 0.3]);
    ylim([1e-11 1e-3]);
    axis square;
    xlabel('电压 (V)');
    ylabel('电流密度 (A)');
    title('电流-电压特性 (对数尺度)');
    legend('Location', 'best');
    grid on;
    

    % Error plot
    figure(errFig);
    nz_idx = data_V ~= 0;
    h_error = bar(data_V(nz_idx), errors(nz_idx), ...
                  'FaceColor', c_data, 'EdgeColor', 'none');
    xlabel('电压 (V)');
    ylabel('相对误差 (%)');
    title(sprintf('拟合误差 (平均: %.2f%%)', avg_error));
    xlim([-0.5 0.3]);
    grid on;
    
    % Parameter display and initial update
    updatePlots(fitFig, errFig, h_fit, h_diode, h_ohmic, h_nonohmic, h_error, ...
        errors, nz_idx, currents, adjusted_params, adjustment_factor, avg_error);
        
    % 持续调整直到用户满意
    while true
        % 显示调整选项
        fprintf('\n当前参数: J0=%.2e, Rs=%.2e, Rsh=%.2e, k=%.2e\n', ...
            adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4));
        fprintf('平均相对误差: %.2f%%\n', avg_error);
        fprintf('\n参数调整选项:\n');
        fprintf('1: 增加 J0  2: 减少 J0\n');
        fprintf('3: 增加 Rs  4: 减少 Rs\n');
        fprintf('5: 增加 Rsh 6: 减少 Rsh\n');
        fprintf('7: 增加 k   8: 减少 k\n');
        fprintf('9: 更改调整步长 (当前: %.2f)\n', adjustment_factor);
        fprintf('0: 结束调整并保存结果\n');
        
        % 获取用户输入并确保是数值类型
        choice_str = input('请选择操作 (0-9): ', 's');
        choice = str2double(choice_str);
        
        % 检查是否为有效数字输入
        if isnan(choice)
            fprintf('请输入有效的数字(0-9)\n');
            continue;
        end
        
        if choice == 0
            currents = calculateCurrents(data_V, adjusted_params, config);
            fit_results.JD = currents.total;
            saveResults(data_V, data_JD, adjusted_params, fit_results, currents);
            saveAdjustedParameters(adjusted_params);
            break;
        elseif choice == 9
            % 调整步长
            new_factor_str = input(sprintf('输入新的调整步长 (当前: %.2f): ', adjustment_factor), 's');
            new_factor = str2double(new_factor_str);
            if ~isnan(new_factor) && new_factor > 0
                adjustment_factor = new_factor;
            else
                fprintf('输入无效，保持当前步长: %.2f\n', adjustment_factor);
            end
            continue;
        elseif choice >= 1 && choice <= 8
              [adjusted_params, currents, errors, avg_error] = applyAdjustment(...
                choice, adjusted_params, adjustment_factor, data_V, data_JD, config, nz_idx);
            updatePlots(fitFig, errFig, h_fit, h_diode, h_ohmic, h_nonohmic, ...
                h_error, errors, nz_idx, currents, adjusted_params, ...
                adjustment_factor, avg_error);
        else
            fprintf('无效的选择，请输入0-9之间的数字\n');
        end
    end
    
    % 最终检查确保Rs为正值
    if adjusted_params(2) <= 0
        fprintf('警告: 最终Rs为负值或零，已调整为正值\n');
        adjusted_params(2) = 1; % 设置为一个合理的小正值
    end
    
    % 计算最终拟合结果
    %fit_results.JD = diodeModel(data_V, adjusted_params, config);
    final_currents = calculateCurrents(data_V, adjusted_params, config);
    fit_results.JD = final_currents.total;
    fit_results.resnorm = sum(((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)).^2);
end
function updatePlots(fitFig, errFig, h_fit, h_diode, h_ohmic, h_nonohmic, h_error, ...
    errors, nz_idx, currents, params, factor, avg_err)
    % Update the plot lines and annotation text.
    set(0,'CurrentFigure',fitFig);
    set(h_fit,      'YData', abs(currents.total));
    set(h_diode,    'YData', abs(currents.diode));
    set(h_ohmic,    'YData', abs(currents.ohmic));
    set(h_nonohmic, 'YData', abs(currents.nonohmic));
    xlim([-0.5 0.3]);
    ylim([1e-11 1e-3]);
    axis square;
    delete(findall(fitFig,'Type','annotation'));
    annotation('textbox',[0.01,0.01,0.98,0.08], ...
        'String',sprintf('J0: %.2e A   Rs: %.2e Ohm   Rsh: %.2e Ohm   k: %.2e   调整步长: %.2f', ...
        params(1),params(2),params(3),params(4),factor), ...
        'EdgeColor','none','FontSize',10,'HorizontalAlignment','center');

    set(0,'CurrentFigure',errFig);
    set(h_error,'YData',errors(nz_idx));
    title(sprintf('拟合误差 (平均: %.2f%%)', avg_err));
    xlim([-0.5 0.3]);
    drawnow;
end

function [params, currents, errors, avg_err] = applyAdjustment(choice, params, factor, data_V, data_JD, config, nz_idx)
    % Apply the user adjustment to parameters and recalc currents/errors.
    param_idx = ceil(choice/2);
    if mod(choice,2)==1
        direction = 1;
    else
        direction = -1;
    end

    delta = params(param_idx) * 0.1 * factor * direction;
    params(param_idx) = params(param_idx) + delta;

    % Enforce reasonable bounds
    if param_idx == 1
        params(param_idx) = max(1e-12, params(param_idx));
    elseif param_idx == 2
        params(param_idx) = max(1, params(param_idx));
        if params(param_idx) <= 0
            fprintf('警告: Rs不能为负值或零。已调整为正值。\n');
            params(param_idx) = 1;
        end
    elseif param_idx == 3
        params(param_idx) = max(1e4, params(param_idx));
    elseif param_idx == 4
        params(param_idx) = max(1e-10, params(param_idx));
    end

    currents = calculateCurrents(data_V, params, config);
    fit_JD = currents.total;
    errors = abs((fit_JD - data_JD) ./ (abs(data_JD)+eps)) * 100;
    avg_err = mean(errors(nz_idx));
end
