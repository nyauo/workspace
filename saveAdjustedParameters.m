function saveAdjustedParameters(params)
    % 生成时间戳
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('adjusted_params_%s.mat', timestamp);
    
    % 保存到文件
    save(filename, 'params');
    fprintf('参数已保存到文件: %s\n', filename);
    
    % 同时导出为文本文件，便于查看
    txt_filename = sprintf('adjusted_params_%s.txt', timestamp);
    fid = fopen(txt_filename, 'w');

    fprintf(fid, 'J01 = %.6e A\n', params(1));
    fprintf(fid, 'Rs = %.6e Ohm\n', params(2));
    fprintf(fid, 'Rsh = %.6e Ohm\n', params(3));
    fprintf(fid, 'k = %.6e\n', params(4));
    fprintf(fid, 'J02 = %.6e A\n', params(5));
    
    fclose(fid);
    fprintf('参数已导出为文本文件: %s\n', txt_filename);
end
