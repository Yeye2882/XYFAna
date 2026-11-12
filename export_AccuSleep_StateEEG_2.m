function fig = export_AccuSleep_StateEEG(dff, neuronFs, neuronStartSec)
% 导出上下时间对齐的 State + EEG + 神经元 DF/F 热图。
% 用法：
% [dff,F0,figDff] = plot_cnmfe_dff_svd(neuron,1.68, ...
%     'OutputPrefix','neuron_dff');
% fig = export_AccuSleep_StateEEG(dff,1.68,0);
% neuronStartSec 是第1帧相对于 EEG 起点的秒数：晚开始为正。
% 不传参数时选择已保存的 neuron_dff.mat（含 dff、Fs、t）。
% 默认假定同步起始；真实同步需根据触发信号确定偏移。
% 整行 NaN/Inf 跳过但保留原编号；局部缺失透明显示，色阶忽略缺失值。
% 本文件须命名为 export_AccuSleep_StateEEG.m。
if nargin < 1 || isempty(dff)
    [fn, fp] = uigetfile('*.mat', '选择已计算的 neuron_dff.mat');
    if isequal(fn,0), error('已取消文件选择。'); end
    D = load(fullfile(fp,fn));
    assert(isfield(D,'dff'), '文件必须含有 dff 矩阵。');
    dff = D.dff;
    if nargin < 2 || isempty(neuronFs)
        assert(isfield(D,'Fs'), '文件没有 Fs，请显式提供神经元采样率。');
        neuronFs = D.Fs;
    end
end
if ~exist('neuronFs','var') || isempty(neuronFs), neuronFs =  2.5463; end
if nargin < 3 || isempty(neuronStartSec), neuronStartSec = 0; end
validateattributes(neuronFs, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(neuronStartSec, {'numeric'}, {'real','finite','scalar'});
% 无效 F0 产生的 NaN 是缺失标记，不能要求整个 dff 都有限。
validateattributes(dff, {'numeric'}, {'real','nonempty','2d'});
dff = double(full(dff));
originalNeuronIDs = 1:size(dff,1);
% 从 MAT 读取时，如保存了 F0，同时屏蔽无效基线对应的整行。
% 不使用 D.ids 重排：此前保存的 dff 包含全部原始行，ids 仅记录绘图选择。
if exist('D','var') && isfield(D,'F0') && ~isempty(D.F0)
    baseline = D.F0;
    assert(isnumeric(baseline) && isreal(baseline) && ismatrix(baseline) && ...
        size(baseline,1)==size(dff,1) && ...
        (size(baseline,2)==1 || size(baseline,2)==size(dff,2)), ...
        '保存的 F0 应为 神经元数 x 1 或 神经元数 x 帧数。');
    invalidBaseline = any(~isfinite(baseline) | baseline<=0,2);
    if any(invalidBaseline)
        warning('export_AccuSleep_StateEEG:InvalidBaseline', ...
            '以下神经元的 F0 无效，将跳过：%s', ...
            mat2str(originalNeuronIDs(invalidBaseline)));
        dff(invalidBaseline,:) = NaN;
    end
end
% Inf 同样视为缺失，不替换成 0，也不对缺失时间点插值。
validNeuron = any(isfinite(dff),2);
skippedNeuronIDs = originalNeuronIDs(~validNeuron);
if ~isempty(skippedNeuronIDs)
    warning('export_AccuSleep_StateEEG:InvalidNeuron', ...
        '以下神经元整行无有效 DF/F，已跳过；其余保留原编号：%s', ...
        mat2str(skippedNeuronIDs));
end
assert(any(validNeuron), '全部神经元均无有效 DF/F，请检查 F0 或 dff。');
neuronIDs = originalNeuronIDs(validNeuron);
dff = dff(validNeuron,:);
missing = ~isfinite(dff);
if any(missing(:))
    warning('export_AccuSleep_StateEEG:MissingSamples', ...
        '保留的神经元中有 %d 个 NaN/Inf 时间点，将显示为空白。', nnz(missing));
    dff(missing) = NaN;
end
assert(size(dff,2)>=2, 'dff 必须是神经元 x 帧矩阵，至少两帧。');
% 与 plot_cnmfe_dff_svd 一致，第1帧 t=0；不插值、不拉伸至 EEG 长度。
tNeuronSec = neuronStartSec + (0:size(dff,2)-1)/neuronFs;
assert(all(isfinite(tNeuronSec)) && all(diff(tNeuronSec)>0), ...
    '神经元时间轴溢出或不递增，请检查 neuronFs 和 neuronStartSec。');

% 原始采样率：1000 Hz；标签 epoch：2.5 s
% 使用 AccuSleep 原有频谱函数，不需要打开 AccuSleep GUI。

%% ================= 用户参数 =================
accuSleepDir = 'F:\Liulab_XYF\Code\Sleep\SleepAnnotation';   % 改成你的 AccuSleep 文件夹

fs       = 1000;
epochSec = 5;

% 导出时间范围，单位：小时
% [] 表示全部；例如 [0.5 0.6] 表示第 30~36 分钟
timeRangeHr = [];

% 输出尺寸与分辨率
widthInch  = 28;
heightInch = 7.2;
dpi        = 600;              % 28 x 7.2 英寸 -> 16800 x 4320 像素；可改 1200

% 频谱亮度：1 为默认；0.8 更亮；1.2 更暗
brightness = 1;

% 若希望自行指定频谱颜色范围，填 [下限 上限]
% 这里使用 AccuSleep 的原始谱值，不是 dB
manualCLim = [];
neuronCLim = [];             % 神经元色阶；[] 使用原脚本的1%~99%规则

%% ================= 准备函数 =================
assert(isfolder(accuSleepDir), '请先修改 accuSleepDir。');
addpath(genpath(accuSleepDir));

assert(exist('createSpectrogram', 'file') == 2, ...
    '未找到 AccuSleep 的 createSpectrogram.m。');
assert(exist('standardizeSR', 'file') == 2, ...
    '未找到 AccuSleep 的 standardizeSR.m。');
assert(exist('AccuSleep_colormap', 'file') == 2, ...
    '未找到 AccuSleep_colormap.m。');
assert(brightness > 0, 'brightness 必须大于 0。');

%% ================= 读取 MAT 文件 =================
[EEG, eegPath] = selectVector('选择 eeg.mat', {'EEG','eeg'});
[EMG, ~]      = selectVector('选择 emg.mat', {'EMG','emg'});
[labels, ~]   = selectVector('选择 label.mat', ...
                             {'labels','label','Labels'});

assert(numel(EEG) == numel(EMG), ...
    'EEG 与 EMG 的采样点数不一致。');
assert(all(isfinite(EEG)) && all(isfinite(EMG)), ...
    'EEG 或 EMG 中含有 NaN/Inf。');
assert(all(ismember(labels, 1:4)), ...
    '标签必须为 1=REM、2=Wake、3=NREM、4=未标注。');

% 本次只绘制 State 和 EEG，EMG 仅用于检查长度
clear EMG

%% ================= 按 AccuSleep 计算频谱 =================
fprintf('正在计算 EEG 频谱...\n');

EEG = EEG - mean(EEG);

% 有意使用原项目的重采样函数，以保持算法一致
eeg128 = standardizeSR(EEG, fs, 128);
clear EEG

[S, ~, F] = createSpectrogram(eeg128, 128, epochSec);
clear eeg128

nEpoch = size(S, 1);

assert(numel(labels) == nEpoch, ...
    ['标签有 %d 个，但频谱有 %d 个 epoch。\n' ...
     '请检查采样率、epoch 长度及标签对应的记录。'], ...
     numel(labels), nEpoch);
assert(nEpoch >= 2, '记录过短。');

F = F(:);
freqKeep = F <= 30;
F = F(freqKeep);
S = S(:, freqKeep);

% 每个 epoch 的中心时刻
timeHr = ((0:nEpoch-1) + 0.5) * epochSec / 3600;
durationHr = nEpoch * epochSec / 3600;

%% ================= 频谱颜色范围 =================
% GUI 随机抽取部分数据计算百分位。
% 这里使用全部频谱，避免每次运行亮度随机变化。
if isempty(manualCLim)
    colorLimits = prctile(S(:), [6 98]);
    if colorLimits(2) <= colorLimits(1)
        colorLimits(2) = colorLimits(1) + ...
                        max(1, abs(colorLimits(1))) * 1e-6;
    end
    colorLimits(2) = colorLimits(1) + ...
        brightness * diff(colorLimits);
else
    assert(numel(manualCLim) == 2 && ...
           all(isfinite(manualCLim)) && ...
           manualCLim(2) > manualCLim(1), ...
           'manualCLim 必须是递增的两个有限数值。');
    colorLimits = manualCLim;
end

%% ================= 导出范围 =================
if isempty(timeRangeHr)
    timeRangeHr = [0 durationHr];
end

assert(numel(timeRangeHr) == 2 && ...
       all(isfinite(timeRangeHr)) && ...
       timeRangeHr(1) >= 0 && ...
       timeRangeHr(2) <= durationHr && ...
       timeRangeHr(2) > timeRangeHr(1), ...
       'timeRangeHr 必须位于 [0, %.6f] 小时内。', durationHr);

%% ================= 创建图窗 =================
fig = figure('Color', 'w', ...
    'Name', 'AccuSleep State + EEG + Neurons', ...
    'NumberTitle', 'off', ...
    'Units', 'inches', ...
    'Position', [0.5 0.5 widthInch heightInch]);

axState = axes('Parent', fig, ...
    'Position', [0.085 0.74 0.90 0.23]);
axEEG = axes('Parent', fig, ...
    'Position', [0.085 0.19 0.90 0.49]);

set([axState axEEG], ...
    'FontName', 'Arial', 'FontSize', 12, ...
    'XColor', 'k', 'YColor', 'k', ...
    'Color', 'w', 'LineWidth', 0.8, ...
    'TickDir', 'out', 'Box', 'off');

%% ================= State：用矢量色块绘制 =================
% 合并连续相同状态，减少 SVG 对象数量
stateColors = [ ...
    0.47  0.67  0.19;        % REM：绿色
    0.14  0.20  0.57;        % Wake：蓝色
    0.996 0.758 0.039;       % NREM：黄色
    0     0     0];          % 未标注：黑色

runStart = [1; find(diff(labels) ~= 0) + 1];
runStop  = [runStart(2:end)-1; nEpoch];

hold(axState, 'on');

for k = 1:numel(runStart)
    state = labels(runStart(k));
    x0 = (runStart(k)-1) * epochSec / 3600;
    x1 = runStop(k) * epochSec / 3600;

    % 只绘制与导出时间段相交的状态块
    x0 = max(x0, timeRangeHr(1));
    x1 = min(x1, timeRangeHr(2));
    if x1 <= x0
        continue
    end

    if state == 4
        y0 = 0.5; y1 = 3.5;
    else
        y0 = state-0.5; y1 = state+0.5;
    end

    patch(axState, [x0 x1 x1 x0], [y0 y0 y1 y1], ...
        stateColors(state,:), 'EdgeColor', 'none');
end

set(axState, 'YLim', [0.5 3.5], 'YDir', 'normal', ...
    'YTick', [1 2 3], 'YTickLabel', {'R','W','N'}, ...
    'XLim', timeRangeHr, 'XTick', []);
ylabel(axState, 'State');

%% ================= EEG 频谱 =================
% 保留完整频谱矩阵；FIG 中可以继续缩放、平移
imagesc(axEEG, timeHr, F, S.');
set(axEEG, 'YDir', 'normal', ...
    'YLim', [0 30], 'YTick', 0:10:30, ...
    'XLim', timeRangeHr, 'CLim', colorLimits);

colormap(axEEG, AccuSleep_colormap());


%% 统一字体和左侧标签布局
set([axState, axEEG], ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...          % 两个坐标轴的刻度字号一致
    'FontWeight', 'normal', ...
    'TickLabelInterpreter', 'none');

% 保留紧凑左边距，右侧为神经元色条留空间，三图绘图区对齐
set(axState, 'Units', 'normalized', ...
    'Position', [0.04 0.84 0.90 0.12]);
set(axEEG, 'Units', 'normalized', ...
    'Position', [0.04 0.58 0.90 0.22]);

% State：保留原来的竖排标签
ylabel(axState, 'State', ...
    'FontName', 'Arial', 'FontSize', 12, ...
    'FontWeight', 'normal', ...
    'Rotation', 90, ...
    'Units', 'normalized', ...
    'Position', [-0.035 0.5 0], ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');

% EEG：保留原来的竖排标签，与 State 对齐
ylabel(axEEG, 'EEG(Hz)', ...
    'FontName', 'Arial', 'FontSize', 12, ...
    'FontWeight', 'normal', ...
    'Rotation', 90, ...
    'Units', 'normalized', ...
    'Position', [-0.035 0.5 0], ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');

% 将频率单位直接写在刻度上，避免标题拥挤
set(axEEG, 'YTick', [0 10 20 30], ...
    'YTickLabel', {'0', '10', '20', '30'});

%% ================= 神经元活动：共同时间坐标 =================
% 三个绘图区采用相同左右边界，colorbar 单独放在右侧。
axNeuron = axes('Parent',fig,'Units','normalized', ...
    'Position',[0.04 0.10 0.90 0.43]);
tNeuronHr = tNeuronSec/3600;
if tNeuronHr(end) < timeRangeHr(1) || tNeuronHr(1) > timeRangeHr(2)
    error('神经元与所选 EEG 时间段无重叠，请检查 neuronStartSec 和 neuronFs。');
end
fprintf('EEG/标签范围：0 ~ %.3f s；神经元帧时间：%.3f ~ %.3f s。\n', ...
    durationHr*3600,tNeuronSec(1),tNeuronSec(end));
if tNeuronHr(1)>timeRangeHr(1) || tNeuronHr(end)<timeRangeHr(2)-1/neuronFs/3600
    warning('神经元未覆盖整个显示时间段；无数据处留白，不拉伸时间。');
end
if isempty(neuronCLim)
    v = sort(dff(isfinite(dff)));
    lo = min(0,v(max(1,ceil(.01*numel(v)))));
    hi = v(max(1,ceil(.99*numel(v))));
    if hi<=lo, hi=lo+1; end
    neuronCLim = [lo hi];
end
validateattributes(neuronCLim,{'numeric'},{'real','finite','numel',2});
assert(neuronCLim(2)>neuronCLim(1),'neuronCLim 必须递增。');
% 显式透明度：缺失点不参与着色；CData 占位值仅用于渲染。
plotData = dff;
finiteMask = isfinite(plotData);
plotData(~finiteMask) = 0;
hNeuron = imagesc(axNeuron,tNeuronHr,1:size(dff,1),plotData);
set(hNeuron,'AlphaData',double(finiteMask),'AlphaDataMapping','none');
% FIG 中保留含 NaN 的数据与原始行号，便于追溯。
setappdata(fig,'NeuronDff',dff);
setappdata(fig,'NeuronIDs',neuronIDs);
setappdata(fig,'SkippedNeuronIDs',skippedNeuronIDs);
colormap(axNeuron,parula(256));
set(axNeuron,'CLim',neuronCLim,'YDir','reverse','Color','w', ...
    'YLim',[.5 size(dff,1)+.5], ...
    'FontName','Arial','FontSize',11,'TickDir','out', ...
    'TickLabelInterpreter','none','Box','off','LineWidth',.8);
step = max(1,ceil(size(dff,1)/8));
tickRows = unique([1 step:step:size(dff,1) size(dff,1)]);
set(axNeuron,'YTick',tickRows,'YTickLabel', ...
    arrayfun(@num2str,neuronIDs(tickRows),'UniformOutput',false));
ylabel(axNeuron,'Neuron ID','FontName','Arial','FontSize',12, ...
    'Units','normalized','Position',[-.035 .5 0]);
cb = colorbar(axNeuron);
set(cb,'Units','normalized','Position',[.95 .10 .009 .43], ...
    'FontName','Arial','FontSize',11);
cb.Label.String = '\Delta F/F_0'; cb.Label.Interpreter = 'tex';
% colorbar 创建后重新固定绘图区位置，保证横向严格对齐。
axNeuron.Position = [0.04 0.10 0.90 0.43];
allAxes = [axState axEEG axNeuron];
if isprop(axState,'PositionConstraint')
    set(allAxes,'PositionConstraint','innerposition');
end
xlabel(axEEG,'');
set(axEEG,'XTickLabel',[]);
xlabel(axNeuron,'Time (hr)','FontName','Arial','FontSize',12);
linkaxes(allAxes,'x');
xlim(axNeuron,timeRangeHr);
drawnow;
% 统一刻度位置；最下方显示时间文字，上方只保留对应刻度。
set(axEEG,'XTick',get(axNeuron,'XTick'));

%% ================= 保存 =================
outDir = fullfile(fileparts(eegPath), 'State_EEG_Neuron_export');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

name = sprintf('State_EEG_Neuron_%.4f-%.4f_hr', ...
               timeRangeHr(1), timeRangeHr(2));
base = fullfile(outDir, name);

set(fig, 'PaperUnits', 'inches', ...
    'PaperPosition', [0 0 widthInch heightInch], ...
    'PaperSize', [widthInch heightInch], ...
    'PaperPositionMode', 'manual', ...
    'InvertHardcopy', 'off');

% FIG 保留数据和图形对象
savefig(fig, [base '.fig']);

% PNG：适合直接插入 PowerPoint
print(fig, [base '.png'], '-dpng', sprintf('-r%d', dpi));

% SVG：State 和文字为矢量，EEG 为嵌入的频谱图像
print(fig, [base '.svg'], '-dsvg', '-painters');

fprintf('\n导出完成，文件夹：\n%s\n', outDir);
fprintf('频谱颜色范围：[%.6g, %.6g]\n', colorLimits);
end

%% ================= 辅助函数 =================
function [v, filePath] = selectVector(dialogTitle, preferredNames)
[fileName, folder] = uigetfile('*.mat', dialogTitle);
if isequal(fileName, 0)
    error('已取消文件选择。');
end

filePath = fullfile(folder, fileName);
data = load(filePath);
names = fieldnames(data);

% 只接受非空、实数、数值向量
valid = cellfun(@(n) isnumeric(data.(n)) && ...
    isreal(data.(n)) && isvector(data.(n)) && ...
    ~isempty(data.(n)), names);
names = names(valid);

assert(~isempty(names), ...
    '%s 中没有数值向量，请检查 MAT 文件结构。', fileName);

selected = '';
for k = 1:numel(preferredNames)
    if any(strcmp(names, preferredNames{k}))
        selected = preferredNames{k};
        break
    end
end

if isempty(selected)
    if numel(names) == 1
        selected = names{1};
    else
        [idx, ok] = listdlg('ListString', names, ...
            'SelectionMode', 'single', ...
            'PromptString', ['选择变量：' fileName]);
        if ~ok
            error('已取消变量选择。');
        end
        selected = names{idx};
    end
end

v = double(data.(selected));
v = v(:);
fprintf('已读取 %s：变量 %s，长度 %d\n', ...
        fileName, selected, numel(v));
end