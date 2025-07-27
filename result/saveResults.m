function saveResults(data_V, data_JD, params, fit_results, currents)
    % 生成时间戳
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    
    outputDir = 'outputs';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    % 保存数据和拟合结果
    results_filename = fullfile(outputDir, sprintf('fit_results_%s.mat', timestamp));
    save(results_filename, 'data_V', 'data_JD', 'params', 'fit_results', 'currents');
    fprintf('拟合结果已保存到文件: %s\n', results_filename);
    
    % 保存图形
    fig_filename = fullfile(outputDir, sprintf('fit_plot_%s.png', timestamp));
    saveas(gcf, fig_filename);
    fprintf('拟合图形已保存到文件: %s\n', fig_filename);
    
    % 导出详细数据到CSV
    csv_filename = fullfile(outputDir, sprintf('fit_data_%s.csv', timestamp));
    fid = fopen(csv_filename, 'w');
    fprintf(fid, 'Voltage(V),Measured_Current(A),Fitted_Current(A),Diode_Current(A),Ohmic_Current(A),Nonohmic_Current(A),Relative_Error(%%)\n');
    
    % 计算相对误差
    rel_errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    
    % 写入数据
    for i = 1:length(data_V)
        fprintf(fid, '%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.2f\n', ...
            data_V(i), data_JD(i), fit_results.JD(i), ...
            currents.diode(i), currents.ohmic(i), currents.nonohmic(i), ...
            rel_errors(i));
    end
    fclose(fid);
    fprintf('拟合数据已导出到CSV文件: %s\n', csv_filename);
    
    % 导出拟合参数到文本文件
    params_filename = fullfile(outputDir, sprintf('fit_params_%s.txt', timestamp));
    fid = fopen(params_filename, 'w');
    fprintf(fid, '拟合参数:\n');
    fprintf(fid, 'J0 = %.6e A\n', params(1));
    fprintf(fid, 'Rs = %.6e Ohm\n', params(2));
    fprintf(fid, 'Rsh = %.6e Ohm\n', params(3));
    fprintf(fid, 'k = %.6e\n\n', params(4));
    
    % 添加误差统计
    fprintf(fid, '拟合误差统计:\n');
    %fprintf(fid, '平均相对误差: %.2f%%\n', mean(rel_errors));
    avg_rel_error = mean(rel_errors(data_V ~= 0));
    fprintf(fid, '平均相对误差: %.2f%%\n', avg_rel_error);
    fprintf(fid, '最大相对误差: %.2f%%\n', max(rel_errors));
    fprintf(fid, '中位相对误差: %.2f%%\n', median(rel_errors));
    
    % 计算不同电压区域的误差
    neg_idx = find(data_V < 0);
    pos_idx = find(data_V > 0);
    if ~isempty(neg_idx)
        fprintf(fid, '负电压区域平均相对误差: %.2f%%\n', mean(rel_errors(neg_idx)));
    end
    if ~isempty(pos_idx)
        fprintf(fid, '正电压区域平均相对误差: %.2f%%\n', mean(rel_errors(pos_idx)));
    end
    
    fclose(fid);
    fprintf('拟合参数和统计信息已保存到文件: %s\n', params_filename);
end
