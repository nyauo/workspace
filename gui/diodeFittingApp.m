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
        ResultsAxes                     matlab.ui.control.UIAxes
        LogTextArea                     matlab.ui.control.TextArea
        NEdit                           matlab.ui.control.NumericEditField
        MEdit                           matlab.ui.control.NumericEditField
        VStartEdit                      matlab.ui.control.NumericEditField
        VStepEdit                       matlab.ui.control.NumericEditField
        VEndEdit                        matlab.ui.control.NumericEditField
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
            app.UIFigure = uifigure('Position',[100 100 450 560], ...
                'Name','Diode Fitting App');
            gl = uigridlayout(app.UIFigure,[12 1]);
            gl.RowHeight = {30,30,30,30,30,40,40,40,40,40,40,'1x'};
            gl.ColumnWidth = {'1x'};
            app.NEdit = uieditfield(gl,'numeric','Value',app.config.physics.n, ...
                'Tooltip','Ideality factor n');
            app.MEdit = uieditfield(gl,'numeric','Value',app.config.physics.m, ...
                'Tooltip','Exponent m');
            app.VStartEdit = uieditfield(gl,'numeric','Value',-0.5, ...
                'Tooltip','Voltage start');
            app.VStepEdit  = uieditfield(gl,'numeric','Value',0.01, ...
                'Tooltip','Voltage step');
            app.VEndEdit   = uieditfield(gl,'numeric','Value',0.3, ...
                'Tooltip','Voltage end');

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

            app.ResultsAxes = uiaxes(gl);
            title(app.ResultsAxes,'Results');
            xlabel(app.ResultsAxes,'Voltage (V)');
            ylabel(app.ResultsAxes,'Current Density (A)');
            app.LogTextArea = uitextarea(gl,'Editable','off');
            guiLog('setHandle', app.LogTextArea);
        end
        function updateConfigFromFields(app)
            app.config.physics.n = app.NEdit.Value;
            app.config.physics.m = app.MEdit.Value;
        end
        function plotOnAxes(app)
            if isempty(app.data_V) || isempty(app.currents)
                cla(app.ResultsAxes);
                return
            end
            ax = app.ResultsAxes;
            cla(ax);
            hold(ax,'on');
            c_data   = [107,174,214]/255; % #6BAED6
            c_total  = [251,106, 74]/255; % #FB6A4A
            c_ohmic  = [144,186, 72]/255; % #90BA48
            c_diode  = [ 19,106,238]/255; % #136AEE
            c_nonohm = [223, 66,227]/255; % #DF42E3
            semilogy(ax, app.data_V, abs(app.data_JD), 'o', 'Color', c_data, 'DisplayName','Measured', 'MarkerSize', 6);
            semilogy(ax, app.data_V, abs(app.currents.total), 'o', 'Color', c_total, 'DisplayName','Total fit', 'MarkerSize', 6);
            semilogy(ax, app.data_V, abs(app.currents.diode), '--', 'Color', c_diode, 'DisplayName','Diode current', 'LineWidth', 1.5);
            semilogy(ax, app.data_V, abs(app.currents.ohmic), '--', 'Color', c_ohmic, 'DisplayName','Ohmic current', 'LineWidth', 1.5);
            semilogy(ax, app.data_V, abs(app.currents.nonohmic), '--', 'Color', c_nonohm, 'DisplayName','Non-ohmic current', 'LineWidth', 1.5);
            xlim(ax,[-0.5 0.3]);
            ylim(ax,[1e-11 1e-3]);
            pbaspect(ax,[1 1 1]);
            grid(ax,'on');
            legend(ax,'Location','best');
            title(ax,'I-V fitting results');
            hold(ax,'off');
        end

        function onLoadData(app)
            try
                app.updateConfigFromFields();
                [app.data_V, app.data_JD] = loadData(app.VStartEdit.Value, ...
                    app.VStepEdit.Value, app.VEndEdit.Value);
            catch ME
                uialert(app.UIFigure, ME.message, 'Data Error');
            end
        end

        function onLoadParams(app)
            outputDir = 'outputs';
            if ~exist(outputDir,'dir')
                mkdir(outputDir);
            end
            [file,path] = uigetfile(fullfile(outputDir,'adjusted_params_*.mat'), ...
                'Choose parameter file');
            if isequal(file,0)
                uialert(app.UIFigure,'No parameter file selected, estimating initial values with Lambert W.','Info');
                if isempty(app.data_V)
                    try
                        [app.data_V, app.data_JD] = loadData(app.VStartEdit.Value, ...
                            app.VStepEdit.Value, app.VEndEdit.Value);
                    catch ME
                        uialert(app.UIFigure, ME.message, 'Data Error');
                        return
                    end
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
                    uialert(app.UIFigure,'File does not contain params variable','Error');
                end
            end
        end

        function onRunFitting(app)
            if isempty(app.data_V)
                try
                    [app.data_V, app.data_JD] = loadData(app.VStartEdit.Value, ...
                        app.VStepEdit.Value, app.VEndEdit.Value);
                catch ME
                    uialert(app.UIFigure, ME.message, 'Data Error');
                    return
                end
            end
            app.updateConfigFromFields();
            validateInputData(app.data_V, app.data_JD);
            if isempty(app.params)
                app.params = initializeParameters(app.data_V, app.data_JD, app.config);
            end
            [app.params, app.fit_results] = performFitting(app.data_V, app.data_JD, app.params, app.config);
            app.currents = calculateCurrents(app.data_V, app.params, app.config);
            plotResults(app.data_V, app.data_JD, app.fit_results, app.currents); % also show external figure
            app.plotOnAxes();
        end

        function onInteractiveAdjust(app)
            if isempty(app.params)
                uialert(app.UIFigure,'Run fitting first.','Info');
                return
            end
            app.updateConfigFromFields();
            [new_params, app.fit_results] = interactiveParameterAdjustment( ...
                app.data_V, app.data_JD, app.params.x0, app.config);
            app.params.x0 = new_params;
            app.currents = calculateCurrents(app.data_V, app.params.x0, app.config);
            app.plotOnAxes();
        end

        function onSaveResults(app)
            if isempty(app.params)
                uialert(app.UIFigure,'No results to save.','Info');
                return
            end
            saveResults(app.data_V, app.data_JD, app.params, app.fit_results, app.currents);
        end

        function onExit(app)
            guiLog('clear');
            delete(app.UIFigure);
        end
    end

    methods (Access = public)
        function app = diodeFittingApp()
            runStartupFcn(app,@startupFcn);
            createComponents(app);
            if nargout == 0
                clear app
            end
        end
    end
end
