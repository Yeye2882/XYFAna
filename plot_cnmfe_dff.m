
[fileName, folderPath] = uigetfile('*.mat', '选择包含 neuron 的 MAT 文件');

if isequal(fileName, 0)
    return;  % 取消选择
end

data = load(fullfile(folderPath, fileName), 'neuron');

if ~isfield(data, 'neuron')
    error('所选 MAT 文件中没有 neuron 变量。');
end

neuron = data.neuron;

[fileName, folderPath] = uigetfile('*.mat', '请选择分块 MAT 文件');

if isequal(fileName, 0)
    disp('已取消选择分块 MAT 文件。');
    return;
end


neuron.P.mat_data = matfile(fullfile(folderPath, fileName));

% 选择与成像记录对应的 EEG 文件（EEG 采样率为 1000 Hz）。
[eegFile, eegFolder] = uigetfile('*.mat', '请选择 EEG MAT 文件', folderPath);
if isequal(eegFile, 0)
    return;
end

eegPath = fullfile(eegFolder, eegFile);
eegVars = who('-file', eegPath);
eegIdx = strcmpi(eegVars, 'EEG');
assert(nnz(eegIdx) == 1, '所选 MAT 文件中应有一个 EEG 或 eeg 变量。');
eegName = eegVars{eegIdx};
eegData = load(eegPath, eegName);
eeg = eegData.(eegName);
validateattributes(eeg, {'numeric'}, {'vector', 'nonempty'});

% frame_range = [首帧, 末帧]，首尾帧都计入。
frameRange = neuron.frame_range;
validateattributes(frameRange, {'numeric'}, {'vector', 'numel', 2, 'finite', 'integer', 'positive'});
assert(frameRange(2) >= frameRange(1), 'neuron.frame_range 的末帧不能小于首帧。');
nFrames = double(frameRange(2)) - double(frameRange(1)) + 1;
Fs = 1000 * nFrames / numel(eeg);  % Hz

[dff,F0,fig] = plot_cnmfe_dff_svd(neuron,Fs,'OutputPrefix','neuron_dff');

