function [dff, labels_dff, Fs, epochSec, t_dff, F0, EEG, EMG, labels] = align_dff_folder(folder)
% 读取一个 session 文件夹，并将睡眠标签映射到每个 dFF 帧。

files = dir(folder);
names = {files(~[files.isdir]).name};
names = names(endsWith(names, '.mat', 'IgnoreCase', true));

data = load(findMat(folder, names, 'neuron_dff'), 'dff', 'F0');
dff = data.dff;
F0 = data.F0;
EEG = loadVector(findMat(folder, names, 'EEG'), 'EEG');
EMG = loadVector(findMat(folder, names, 'EMG'), 'EMG');
labels = loadVector(findMat(folder, names, 'labels'), 'labels');

% 根据 EEG 时长选择最接近的标签时长，并计算 dFF 采样率
FsEEG = 1000;
nEEG = numel(EEG);
nFrames = size(dff, 2);
assert(nFrames > 0, 'dff 没有可用的帧。');
epochOptions = [2.5, 5];
[~, k] = min(abs(nEEG / FsEEG / numel(labels) - epochOptions));
epochSec = epochOptions(k);
Fs = FsEEG * nFrames / nEEG;

% 第一帧对应 t = 0，每个标签覆盖一个左闭右开的时间区间
frame = 0:nFrames-1;
t_dff = frame / Fs;
idx = floor((frame * nEEG) / (nFrames * FsEEG * epochSec)) + 1;
valid = idx <= numel(labels);
labels_dff = nan(1, nFrames);  % 超出标签覆盖范围的帧保留为 NaN
labels_dff(valid) = labels(idx(valid));

end

function path = findMat(folder, names, key)
    exact = strcmpi(names, [key '.mat']);
    if any(exact)
        match = exact;
    else
        match = contains(names, key, 'IgnoreCase', true);
    end
    assert(nnz(match) == 1, '关键词 "%s" 应匹配一个 MAT 文件，实际匹配到 %d 个。', key, nnz(match));
    path = fullfile(folder, names{match});
end

function v = loadVector(path, key)
    s = load(path);
    fields = fieldnames(s);
    match = strcmpi(fields, key);
    if strcmp(key, 'labels') && ~any(match)
        match = strcmpi(fields, 'label');
    end
    if ~any(match)
        match = cellfun(@(f) isDataVector(s.(f)), fields);
    end
    assert(nnz(match) == 1, '无法在 %s 中确定唯一的 %s 向量。', path, key);
    v = s.(fields{match});
    assert(isDataVector(v), '%s 中的 %s 必须为非空数值或逻辑向量。', path, key);
    v = v(:);
end

function tf = isDataVector(v)
    tf = (isnumeric(v) || islogical(v)) && isvector(v) && ~isempty(v);
end
