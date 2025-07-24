function plotResults(V, JD_measured, fit_results, currents)
    % Figure for semilog I-V curve
    fitFig = figure('Name','拟合结果','Position',[100 100 600 600]);

    % 绘制 I-V 特性（对数坐标）及各电流分量
    figure(fitFig);
    
    % 定义颜色
    c_data   = [107,174,214]/255; % #6BAED6
    c_total  = [251,106, 74]/255; % #FB6A4A
    c_ohmic  = [144,186, 72]/255; % #90BA48
    c_diode  = [ 19,106,238]/255; % #136AEE
    c_tunnel = [172,209,230]/255; % #ACD1E6
    c_nonohm = [223, 66,227]/255; % #DF42E3

    % 绘制测量数据和拟合结果
    semilogy(V, abs(JD_measured), 'o', 'Color', c_data, 'DisplayName', '测量数据', 'MarkerSize', 6);
    hold on;
    semilogy(V, abs(currents.total), 'o', 'Color', c_total, 'DisplayName', '总拟合电流', 'MarkerSize', 6);
    semilogy(V, abs(currents.diode), '--', 'Color', c_diode, 'DisplayName', '二极管电流', 'LineWidth', 1.5);
    semilogy(V, abs(currents.tunnel), '--', 'Color', c_tunnel, 'DisplayName', '隧穿电流', 'LineWidth', 1.5);
    semilogy(V, abs(currents.ohmic), '--', 'Color', c_ohmic, 'DisplayName', '欧姆电流', 'LineWidth', 1.5);
    semilogy(V, abs(currents.nonohmic), '--', 'Color', c_nonohm, 'DisplayName', '非欧姆电流', 'LineWidth', 1.5);
    xlim([-0.5 0.3]);
    ylim([1e-11 1e-3]);
    axis square;
    
    % 设置坐标轴
    xlabel('电压 (V)', 'FontSize', 12);
    ylabel('电流密度 (A)', 'FontSize', 12);
    title('二极管I-V特性曲线及电流分量（对数坐标）', 'FontSize', 14);
    grid on;
    legend('Location', 'best');
    
    % 第二个图形：误差和线性 I-V
    errFig = figure('Name','误差与线性I-V','Position',[750 100 600 600]);
    subplot(2,1,1);
    
    relative_error = abs((currents.total - JD_measured) ./ (abs(JD_measured) + eps)) * 100;
    nz_idx = V ~= 0;  % 忽略零电压点进行误差计算和绘图
       bar(V(nz_idx), relative_error(nz_idx), 'FaceColor', c_data, 'EdgeColor', 'none');
    avg_rel_err = mean(relative_error(nz_idx));
    max_rel_err = max(relative_error(nz_idx));
    xlabel('电压 (V)');
    ylabel('相对误差 (%)');
    title('拟合相对误差');
    xlim([-0.5 0.3]);
    grid on;
    
    % 第三个子图：I-V特性（线性坐标）
    subplot(2,1,2);
    plot(V, JD_measured, 'Color', c_data, 'LineWidth', 1.5, 'DisplayName', '测量数据');
    hold on;
    plot(V, currents.total, 'Color', c_total, 'LineWidth', 1.5, 'DisplayName', '拟合结果');
    xlabel('电压 (V)');
    ylabel('电流密度 (A)');
    title('I-V特性曲线（线性坐标）');
    grid on;
    legend('Location', 'best');
    
    figure(errFig);
    sgtitle('二极管特性拟合结果分析','FontSize',14);

    annotation(errFig,'textbox',[0.15,0.02,0.7,0.05], ...
        'String', sprintf('最大相对误差: %.2f%%  平均相对误差: %.2f%%', max_rel_err, avg_rel_err), ...
        'EdgeColor','none','HorizontalAlignment','center');
    fprintf('二极管电流占比: %.2f%%\n', mean(abs(currents.diode ./ currents.total)) * 100);
    fprintf('欧姆电流占比: %.2f%%\n', mean(abs(currents.ohmic ./ currents.total)) * 100);
    fprintf('隧穿电流占比: %.2f%%\n', mean(abs(currents.tunnel ./ currents.total)) * 100);
    fprintf('非欧姆电流占比: %.2f%%\n', mean(abs(currents.nonohmic ./ currents.total)) * 100);
end
