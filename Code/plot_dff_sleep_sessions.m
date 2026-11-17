function [sessionMeans, fig, mouseMeans, mouseFig] = plot_dff_sleep_sessions()
% 在同一窗口多选“小鼠名_日期_N”文件夹，按回车确认。
% sessionMeans 和 fig 按 session 汇总，mouseMeans 和 mouseFig 按小鼠汇总。

stateNames = {'Wake', 'NREM', 'REM'};
stateCodes = [2 3 1];  % AccuSleep 默认编码
sessionMeans = table('Size', [0 9], ...
    'VariableTypes', {'string','string','string','double','double','double','double','double','double'}, ...
    'VariableNames', {'Session','Mouse','Date','N','Fs','EpochSec','Wake','NREM','REM'});
fig = [];
mouseFig = [];
mouseMeans = table(strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', [{'Mouse', 'NumSessions'}, stateNames]);
sessionPattern = '^(.+)_(\d{8})_(\d+)$';
folders = selectSessionFolders(pwd, sessionPattern);

for k = 1:numel(folders)
    folder = folders{k};
    [~, name] = fileparts(folder);
    parts = regexp(name, sessionPattern, 'tokens', 'once');
    sessionNumber = str2double(parts{3});
    alreadyAdded = sessionMeans.Mouse == string(parts{1}) & ...
        sessionMeans.Date == string(parts{2}) & sessionMeans.N == sessionNumber;
    if any(alreadyAdded)
        warning('Session %s 已添加。', name);
        continue;
    end

    try
        [dff, labels_dff, Fs, epochSec] = align_dff_folder(folder);
    catch err
        warning('无法处理 %s。%s', name, err.message);
        continue;
    end

    % 先对每个 ROI 的状态内帧取均值，再对 ROI 等权平均。
    values = nan(1, 3);
    for s = 1:3
        mask = labels_dff == stateCodes(s);
        if ~any(mask), continue; end
        frames = double(dff(:, mask));
        frames(~isfinite(frames)) = NaN;
        roiMean = mean(frames, 2, 'omitnan');
        values(s) = mean(roiMean, 'omitnan');
    end
    if all(isnan(values))
        warning('Session %s 在三个状态下均无有效 dFF，已跳过。', name);
        continue;
    end
    sessionMeans(end+1, :) = {string(name), string(parts{1}), string(parts{2}), ...
        sessionNumber, Fs, epochSec, values(1), values(2), values(3)};
end

if isempty(sessionMeans), return; end

sessionValues = sessionMeans{:, stateNames};
fig = plotStateMeans(sessionValues, sessionMeans.Session, stateNames, ...
    'Calcium activity across sessions');

% 以文件夹名中的小鼠名分组，每只小鼠内的 session 等权。
[mice, ~, group] = unique(sessionMeans.Mouse, 'stable');
mouseValues = nan(numel(mice), numel(stateNames));
numSessions = accumarray(group, 1);
for m = 1:numel(mice)
    mouseValues(m, :) = mean(sessionValues(group == m, :), 1, 'omitnan');
end
mouseMeans = [table(mice, numSessions, 'VariableNames', {'Mouse', 'NumSessions'}), ...
    array2table(mouseValues, 'VariableNames', stateNames)];
mouseFig = plotStateMeans(mouseValues, mice, stateNames, 'Calcium activity across mice');
end

function fig = plotStateMeans(values, names, stateNames, figureName)
% 每行等权。浅灰线为单个 session 或小鼠，黑线和阴影为均值 ± SEM。
n = sum(isfinite(values), 1);
center = mean(values, 1, 'omitnan');
sem = std(values, 0, 1, 'omitnan') ./ sqrt(n);
sem(n < 2) = NaN;

fig = figure('Color', 'w', 'Name', figureName, ...
    'Position', [200 200 360 400]);
ax = axes('Parent', fig);
hold(ax, 'on');
for k = 1:size(values, 1)
    plot(ax, 1:3, values(k, :), '.-', 'Color', [0.73 0.73 0.73], ...
        'LineWidth', 0.8, 'MarkerSize', 7, 'DisplayName', char(names(k)));
end

hasBand = false(1, 3);
for s = 1:2
    j = s:s+1;
    if all(isfinite(sem(j)))
        fill(ax, [j fliplr(j)], [center(j)-sem(j), fliplr(center(j)+sem(j))], ...
            [0.28 0.28 0.28], 'FaceAlpha', 0.42, 'EdgeColor', 'none');
        hasBand(j) = true;
    end
end
unbanded = isfinite(sem) & ~hasBand;
if any(unbanded)
    errorbar(ax, find(unbanded), center(unbanded), sem(unbanded), 'k', 'LineStyle', 'none');
end
plot(ax, 1:3, center, 'k.-', 'LineWidth', 1.65, 'MarkerSize', 12);
set(ax, 'XTick', 1:3, 'XTickLabel', stateNames, 'TickDir', 'out', ...
    'Box', 'off', 'FontSize', 11, 'LineWidth', 0.9);
xlim(ax, [0.8 3.2]);
ylabel(ax, '\DeltaF/F_0');
end

function folders = selectSessionFolders(currentFolder, sessionPattern)
% 用多选列表浏览同一父目录下的 session。
folders = {};
names = {};
win = dialog('Name', '选择 session 文件夹', 'Position', [100 100 640 430], ...
    'WindowKeyPressFcn', @onKey);
movegui(win, 'center');
pathBox = uicontrol(win, 'Style', 'edit', 'Position', [16 382 460 28], ...
    'HorizontalAlignment', 'left', 'FontSize', 10);
uicontrol(win, 'Style', 'pushbutton', 'String', '转到', ...
    'Position', [486 382 62 28], 'Callback', @(~,~) showFolder(strtrim(pathBox.String)));
uicontrol(win, 'Style', 'pushbutton', 'String', '上一级', ...
    'Position', [558 382 66 28], 'Callback', @(~,~) showFolder(fileparts(currentFolder)));
list = uicontrol(win, 'Style', 'listbox', 'Min', 0, 'Max', 2, ...
    'Position', [16 89 608 282], 'FontSize', 11, 'Callback', @onList);
status = uicontrol(win, 'Style', 'text', 'Position', [16 53 608 25], ...
    'HorizontalAlignment', 'left', 'FontSize', 10);
uicontrol(win, 'Style', 'pushbutton', 'String', '确认', ...
    'Position', [448 15 80 30], 'Callback', @confirm);
uicontrol(win, 'Style', 'pushbutton', 'String', '取消', ...
    'Position', [544 15 80 30], 'Callback', @(~,~) delete(win));

showFolder(currentFolder);
uiwait(win);

    function showFolder(path)
        if ~isfolder(path)
            status.String = '目录不存在，请输入有效的父目录路径。';
            return;
        end
        [ok, info] = fileattrib(path);
        if ~ok
            status.String = '无法读取该目录。';
            return;
        end
        currentFolder = info.Name;
        entries = dir(currentFolder);
        names = sort(setdiff({entries([entries.isdir]).name}, {'.','..'}));
        set(list, 'Value', [], 'String', names);
        pathBox.String = currentFolder;
        status.String = 'Ctrl/Shift 多选 session，回车确认。双击可进入子目录。';
        uicontrol(list);
    end

    function onList(~, ~)
        if strcmp(win.SelectionType, 'open') && isscalar(list.Value)
            showFolder(fullfile(currentFolder, names{list.Value}));
        end
    end

    function onKey(~, event)
        switch event.Key
            case {'return', 'enter'}
                confirm();
            case 'escape'
                delete(win);
        end
    end

    function confirm(varargin)
        chosen = names(list.Value);
        if isempty(chosen)
            status.String = '请至少选择一个 session 文件夹。';
            return;
        end
        if any(cellfun('isempty', regexp(chosen, sessionPattern, 'once')))
            status.String = '请选择名称为“小鼠名_日期_N”的文件夹。';
            return;
        end
        folders = fullfile(currentFolder, chosen);
        delete(win);
    end
end
