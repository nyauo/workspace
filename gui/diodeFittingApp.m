classdef diodeFittingApp < matlab.apps.AppBase
    % Simple GUI wrapper around main fitting workflow

    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        LoadDataButton                  matlab.ui.control.Button
        LoadParamsButton                matlab.ui.control.Button
        RunFittingButton                matlab.ui.control.Button
        InteractiveAdjustButton         matlab.ui.control.Button
        SaveResultsButton               matlab.ui.control.Button
        ExitButton                      matlab.ui.control.Button
    end

    properties (Access = private)
        data_V
        data_JD
        params
        fit_results
        currents
        config
    end

    methods (Access = private)
        function startupFcn(app)
            % Add paths for IO, fitting and result utilities
            appDir = fileparts(mfilename('fullpath'));
            repoDir = fileparts(appDir);
            addpath(genpath(fullfile(repoDir,'io')));
            addpath(genpath(fullfile(repoDir,'fit')));
            addpath(genpath(fullfile(repoDir,'result')));

            app.config = loadConfig();
        end

        function createComponents(app)
            app.UIFigure = uifigure('Position',[100 100 250 260], ...
                'Name','Diode Fitting App');
            gl = uigridlayout(app.UIFigure,[6 1]);
            gl.RowHeight = repmat({'1x'},1,6);
            gl.ColumnWidth = {'1x'};

            app.LoadDataButton = uibutton(gl,'Text','Load Data', ...
                'ButtonPushedFcn',@(src,event)app.onLoadData());
            app.LoadParamsButton = uibutton(gl,'Text','Load Parameters', ...
                'ButtonPushedFcn',@(src,event)app.onLoadParams());
            app.RunFittingButton = uibutton(gl,'Text','Run Fitting', ...
                'ButtonPushedFcn',@(src,event)app.onRunFitting());
            app.InteractiveAdjustButton = uibutton(gl,'Text','Interactive Adjust', ...
                'ButtonPushedFcn',@(src,event)app.onInteractiveAdjust());
            app.SaveResultsButton = uibutton(gl,'Text','Save Results', ...
                'ButtonPushedFcn',@(src,event)app.onSaveResults());
            app.ExitButton = uibutton(gl,'Text','Exit', ...
                'ButtonPushedFcn',@(src,event)app.onExit());
        end

        function onLoadData(app)
            [app.data_V, app.data_JD] = loadData();
        end

        function onLoadParams(app)
            outputDir = 'outputs';
            if ~exist(outputDir,'dir')
                mkdir(outputDir);
            end
            [file,path] = uigetfile(fullfile(outputDir,'adjusted_params_*.mat'), ...
                '选择历史参数文件');
            if isequal(file,0)
                uialert(app.UIFigure,'未选择参数文件，将使用Lambert W估计初始值','提示');
                if isempty(app.data_V)
                    [app.data_V, app.data_JD] = loadData();
                end
                app.params = initializeParameters(app.data_V, app.data_JD, app.config);
            else
                S = load(fullfile(path,file));
                if isfield(S,'params')
                    app.params.x0 = S.params;
                    app.params.ub = [1e-6,1e4,1e10,1e-5];
                    app.params.lb = [1e-12,1e1,1e5,1e-10];
                    app.params.scaleFactors = [1e-9,1e3,1e7,1e-8];
                else
                    uialert(app.UIFigure,'文件不包含 params 变量','错误');
                end
            end
        end

        function onRunFitting(app)
            if isempty(app.data_V)
                [app.data_V, app.data_JD] = loadData();
            end
            validateInputData(app.data_V, app.data_JD);
            if isempty(app.params)
                app.params = initializeParameters(app.data_V, app.data_JD, app.config);
            end
            [app.params, app.fit_results] = performFitting(app.data_V, app.data_JD, app.params, app.config);
            app.currents = calculateCurrents(app.data_V, app.params, app.config);
            plotResults(app.data_V, app.data_JD, app.fit_results, app.currents);
        end

        function onInteractiveAdjust(app)
            if isempty(app.params)
                uialert(app.UIFigure,'请先运行拟合','提示');
                return
            end
            [app.params, app.fit_results] = interactiveParameterAdjustment( ...
                app.data_V, app.data_JD, app.params, app.config);
            app.currents = calculateCurrents(app.data_V, app.params, app.config);
        end

        function onSaveResults(app)
            if isempty(app.params)
                uialert(app.UIFigure,'暂无结果可保存','提示');
                return
            end
            saveResults(app.data_V, app.data_JD, app.params, app.fit_results, app.currents);
        end

        function onExit(app)
            delete(app.UIFigure);
        end
    end

    methods (Access = public)
        function app = diodeFittingApp()
            createComponents(app);
            runStartupFcn(app,@startupFcn);
            if nargout == 0
                clear app
            end
        end
    end
end
