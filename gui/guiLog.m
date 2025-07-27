function guiLog(action, arg)
%GUILOG  Simple GUI logger. Append messages to a text area if available.
%   guiLog(MSG) appends MSG to the registered text area or prints to stdout.
%   guiLog('setHandle', H) registers uitextarea H for logging.
%   guiLog('clear') clears the registered handle.

    persistent logHandle

    if nargin == 0
        return;
    end

    if ischar(action) && (strcmp(action,'setHandle') || strcmp(action,'set'))
        logHandle = arg;
        if isempty(logHandle)
            return;
        end
        if isvalid(logHandle) && isempty(logHandle.Value)
            logHandle.Value = {};
        end
        return;
    elseif ischar(action) && strcmp(action,'clear')
        logHandle = [];
        return;
    end

    msg = action;
    if ~isempty(logHandle) && isvalid(logHandle)
        logHandle.Value = [logHandle.Value; {char(msg)}];
        drawnow;
    else
        fprintf('%s\n', msg);
    end
end
