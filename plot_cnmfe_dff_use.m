
[fileName, folderPath] = uigetfile('*.mat', '选择包含 neuron 的 MAT 文件');

if isequal(fileName, 0)
    return;  % 取消选择
end

data = load(fullfile(folderPath, fileName), 'neuron');

if ~isfield(data, 'neuron')
    error('所选 MAT 文件中没有 neuron 变量。');
end

neuron = data.neuron;

[fileName, folderPath] = uigetfile('*.mat', '请选择数据 MAT 文件');

if isequal(fileName, 0)
    disp('已取消选择 MAT 文件。');
    return;
end

neuron.P.mat_data = matfile(fullfile(folderPath, fileName));

%[dff,F0,fig] = plot_cnmfe_dff_svd(neuron,1.68,'OutputPrefix','neuron_dff');

[dff,F0,fig] = plot_cnmfe_dff_svd(neuron,2.5415,'OutputPrefix','neuron_dff');